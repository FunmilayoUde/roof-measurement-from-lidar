"""Stage 3: facet polygons and typed edges for one house.

    uv run python scripts/stage3_topology.py --house-id delft_03

Reads data/interim/<house-id>/{roof.ply, stats.json}. Segments with the stage 2 choice
(RANSAC, threshold 0.12 m, min_area 1 m²), then finds internal edges and facet polygons.
Writes data/interim/<house-id>/stage3/{topology.json, plan.png}.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import open3d as o3d

from roofrecon.planes import azimuth_deg, hull_area, pitch_deg
from roofrecon.polygons import roof_polygons
from roofrecon.segment import segment_ransac
from roofrecon.topology import internal_edges

# Validated categorical slots (reference palette, light mode), one per edge kind, fixed order.
KIND_COLORS = {"ridge": "#2a78d6", "hip": "#eb6834", "valley": "#1baf7a", "eave": "#eda100",
               "rake": "#e87ba4", "break": "#008300", "top": "#6250d6", "perimeter": "#e34948",
               "other": "#52514e"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--threshold", type=float, default=0.12)
    ap.add_argument("--min-area", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    house = args.root / args.house_id
    xyz = np.asarray(o3d.io.read_point_cloud(str(house / "roof.ply")).points)
    density = json.loads((house / "stats.json").read_text())["density_inside_walls"]["nominal_pts_per_m2"]
    seg = segment_ransac(xyz, threshold=args.threshold, min_points=max(3, round(args.min_area * density)), seed=args.seed)
    edges = internal_edges(xyz, seg)
    polys = roof_polygons(xyz, seg, edges)

    lengths = defaultdict(float)
    facets = []
    for P in polys:
        pts = xyz[seg.labels == P.facet]
        hull = hull_area(pts, P.plane)
        facets.append({
            "facet": P.facet, "points": int(len(pts)), "pitch_deg": round(pitch_deg(P.plane.normal), 2),
            "azimuth_deg": None if (az := azimuth_deg(P.plane.normal)) is None else round(az, 1),
            "polygon_area_m2": round(P.area, 2), "hull_area_m2": round(hull, 2), "plan_area_m2": round(P.plan_area, 2),
            "fallback": P.fallback,
            "sides": [{"kind": s.kind, "neighbour": s.neighbour, "length_m": round(s.length, 2)} for s in P.sides],
            "vertices_local": np.round(P.vertices, 3).tolist(),
        })
        for s in P.sides:
            lengths[s.kind] += s.length / (2 if s.neighbour is not None else 1)

    out = house / "stage3"
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "house_id": args.house_id,
        "planes": len(seg.planes),
        "coverage": round(seg.coverage, 4),
        "fallbacks": sum(P.fallback for P in polys),
        "roof_area_polygon_m2": round(sum(P.area for P in polys), 1),
        "roof_area_hull_m2": round(sum(f["hull_area_m2"] for f in facets), 1),
        "lengths_m": {k: round(v, 2) for k, v in sorted(lengths.items())},
        "internal_edges": [{"kind": e.kind, "facets": list(e.facets), "length_m": round(e.length, 2),
                            "snapped": list(e.snapped)} for e in edges],
        "facets": facets,
    }
    (out / "topology.json").write_text(json.dumps(summary, indent=2))
    plot_plan(xyz, seg.labels, polys, out / "plan.png", f"{args.house_id}: {len(polys)} facets")

    print(f"{args.house_id}: {len(polys)} facets, coverage {seg.coverage:.3f}, fallbacks {summary['fallbacks']}")
    print(f"  {'facet':>5} {'pitch°':>7} {'az°':>6} {'polygon m²':>11} {'hull m²':>8}  sides")
    for f in sorted(facets, key=lambda f: -f["polygon_area_m2"]):
        az = "-" if f["azimuth_deg"] is None else f"{f['azimuth_deg']:.0f}"
        sides = ", ".join(f"{s['kind']} {s['length_m']}" for s in f["sides"])
        print(f"  {f['facet']:>5} {f['pitch_deg']:>7.1f} {az:>6} {f['polygon_area_m2']:>11.1f} {f['hull_area_m2']:>8.1f}  {sides}")
    print(f"  roof area: polygons {summary['roof_area_polygon_m2']} m², hulls {summary['roof_area_hull_m2']} m²")
    print("  lengths: " + ", ".join(f"{k} {v}" for k, v in summary["lengths_m"].items()))


def plot_plan(xyz, labels, polys, path: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from roofrecon.viz import segment_colors

    fig, ax = plt.subplots(figsize=(9, 8), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.scatter(xyz[:, 0], xyz[:, 1], s=2, c=segment_colors(labels), alpha=0.35)
    used = set()
    for P in polys:
        for s in P.sides:
            ax.plot([s.p0[0], s.p1[0]], [s.p0[1], s.p1[1]], color=KIND_COLORS.get(s.kind, "#52514e"), lw=2.2,
                    ls="--" if P.fallback else "-")
            used.add(s.kind)
        c = P.vertices.mean(axis=0)
        ax.text(c[0], c[1], f"{P.facet}\n{P.area:.1f} m²", fontsize=8, ha="center", va="center", color="#0b0b0b",
                bbox={"facecolor": "#fcfcfb", "alpha": 0.8, "pad": 1, "edgecolor": "none"})
    ax.set_aspect("equal")
    ax.set_title(title, loc="left", color="#0b0b0b")
    ax.tick_params(colors="#52514e", labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    handles = [Line2D([], [], color=KIND_COLORS[k], lw=2.2, label=k) for k in KIND_COLORS if k in used]
    if any(P.fallback for P in polys):
        handles.append(Line2D([], [], color="#52514e", lw=2.2, ls="--", label="fallback outline"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=len(handles), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor="#fcfcfb", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
