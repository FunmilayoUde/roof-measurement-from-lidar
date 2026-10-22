"""Stage 4: write the measurement report for one house.

    uv run python scripts/stage4_report.py --house-id delft_03

Reads data/interim/<house-id>/{roof.ply, stats.json, footprint.geojson}, builds the roof
model with the chosen stage 2/3 settings, and writes data/interim/<house-id>/report.json.
Coordinates in the report are local (stage 1 recentering); add source.offset_xyz to get
EPSG:28992 / NAP.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import open3d as o3d

sys.path.insert(0, str(Path(__file__).parent))
from roofrecon.measure import measurement_report  # noqa: E402
from roofrecon.pipeline import build_roof  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    house = args.root / args.house_id
    xyz = np.asarray(o3d.io.read_point_cloud(str(house / "roof.ply")).points)
    stats = json.loads((house / "stats.json").read_text())
    footprint = json.loads((house / "footprint.geojson").read_text())
    meta = {
        "house_id": args.house_id,
        "bag_id": footprint.get("properties", {}).get("identificatie"),
        "tile": Path(stats["tile"]["path"]).name,
        "crs": stats["tile_crs"],
        "offset_xyz": stats["offset_xyz"],
        "footprint_area_m2": round(stats["footprint_area_m2"], 1),
    }
    rep = measurement_report(build_roof(xyz, stats["density_inside_walls"]["nominal_pts_per_m2"], seed=args.seed), meta)
    (house / "report.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False))
    plot_report_plan(rep, house / "report_plan.png")

    t = rep["totals"]
    print(
        f"{args.house_id}: {t['facets']} facets, {t['roof_area_m2']} ± {t['roof_area_p90_error_m2']} m² (p90) "
        f"= {t['roof_area_ft2']:.0f} ft² = {t['squares']} squares"
    )
    print(f"  {'facet':>5} {'pitch':>12} {'facing':>6} {'area m²':>8} {'±p90':>6}  confidence  outline")
    for f in sorted(rep["facets"], key=lambda f: -f["area_m2"]):
        print(
            f"  {f['id']:>5} {f['pitch_deg']:>5}° ({f['pitch_rise_per_12']:>4}/12) "
            f"{f['facing'] or 'flat':>6} {f['area_m2']:>8} "
            f"{f['area_error_p90_pct']:>5}%  {f['confidence']:10}  {f['outline']}"
        )
    print(
        "  lengths: "
        + ", ".join(f"{k} {v} m ({t['length_ft_by_type'][k]} ft)" for k, v in t["length_m_by_type"].items())
    )
    for note in rep["quality"]["notes"]:
        print(f"  note: {note}")


def plot_report_plan(rep: dict, path: Path) -> None:
    """Plan of the report with every edge and facet numbered, for matching hand measurements
    (validation/PROTOCOL.md). Look at it only after measuring."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from stage3_topology import KIND_COLORS

    colors = KIND_COLORS | {"unclassified": "#52514e"}
    fig, ax = plt.subplots(figsize=(10, 9), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    used = set()
    for f in rep["facets"]:
        v = np.asarray(f["vertices_local_m"])
        c = v.mean(axis=0)
        ax.text(c[0], c[1], f"F{f['id']}\n{f['pitch_deg']}°", fontsize=9, ha="center", va="center", color="#52514e")
    for e in rep["edges"]:
        a, b = np.asarray(e["p0_local_m"]), np.asarray(e["p1_local_m"])
        ax.plot([a[0], b[0]], [a[1], b[1]], color=colors.get(e["type"], "#52514e"), lw=2)
        m = (a + b) / 2
        ax.text(
            m[0],
            m[1],
            str(e["id"]),
            fontsize=8,
            ha="center",
            va="center",
            color="#0b0b0b",
            bbox={"facecolor": "#fcfcfb", "alpha": 0.85, "pad": 0.8, "edgecolor": "none"},
        )
        used.add(e["type"])
    ax.set_aspect("equal")
    ax.set_title(
        f"{rep['house_id']}: edge ids (numbers) and facet ids (F) for matching hand measurements",
        loc="left",
        fontsize=10,
        color="#0b0b0b",
    )
    ax.tick_params(colors="#52514e", labelsize=8)
    handles = [Line2D([], [], color=colors[k], lw=2, label=k) for k in colors if k in used]
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=len(handles), frameon=False, fontsize=8
    )
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor="#fcfcfb", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
