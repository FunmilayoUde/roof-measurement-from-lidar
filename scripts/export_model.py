"""Export one house's 3D roof model from its report.

    uv run python scripts/export_model.py --house-id delft_03

Reads data/interim/<house-id>/report.json; writes model.obj (+ .mtl), model.ply, edges.ply and
model_preview.png in data/interim/<house-id>/model/. View with scripts/view_model.py, or open
model.obj in CloudCompare or Blender (local coordinates; the offset to EPSG:28992 is in the
OBJ header).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from roofrecon.export import EDGE_COLORS, FACET_COLORS, mesh_area, roof_mesh, write_obj, write_ply


def preview(report: dict, mesh: dict, xyz: np.ndarray, path: Path) -> None:
    """Static 3D view: facets coloured, edges by type, a thin sample of the laser points."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    ids = [f["id"] for f in report["facets"] if f["id"] in set(mesh["facet_of_face"].tolist())]
    colour_of = {fid: FACET_COLORS[i % len(FACET_COLORS)] for i, fid in enumerate(ids)}
    fig = plt.figure(figsize=(14, 7), facecolor="#fcfcfb")
    for k, (elev, azim) in enumerate([(35, -60), (20, 30)]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        ax.set_facecolor("#fcfcfb")
        tris = mesh["vertices"][mesh["faces"]]
        ax.add_collection3d(
            Poly3DCollection(
                tris, facecolors=[colour_of[f] for f in mesh["facet_of_face"]], edgecolors="none", alpha=0.85
            )
        )
        used = set()
        for e in report["edges"]:
            a, b = np.asarray(e["p0_local_m"]), np.asarray(e["p1_local_m"])
            ax.plot(*zip(a, b), color=EDGE_COLORS.get(e["type"], "#52514e"), lw=2.2)
            used.add(e["type"])
        sample = xyz[:: max(1, len(xyz) // 1500)]
        ax.scatter(*sample.T, s=0.5, c="#52514e", alpha=0.35)
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        span = (hi - lo).max() / 2
        mid = (hi + lo) / 2
        ax.set_xlim(mid[0] - span, mid[0] + span)
        ax.set_ylim(mid[1] - span, mid[1] + span)
        ax.set_zlim(lo[2] - 0.5, lo[2] - 0.5 + 2 * span)
        ax.view_init(elev=elev, azim=azim)
        ax.set_box_aspect((1, 1, 1))
        ax.tick_params(labelsize=7, colors="#52514e")
    handles = [Line2D([], [], color=EDGE_COLORS[t], lw=2.2, label=t) for t in EDGE_COLORS if t in used]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9)
    t = report["totals"]
    fig.suptitle(
        f"{report['house_id']}: {t['facets']} facets, {t['roof_area_m2']} m², laser points in grey",
        x=0.02,
        ha="left",
        color="#0b0b0b",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(path, dpi=110, facecolor="#fcfcfb")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    import open3d as o3d

    house = args.root / args.house_id
    report = json.loads((house / "report.json").read_text())
    xyz = np.asarray(o3d.io.read_point_cloud(str(house / "roof.ply")).points)
    out = house / "model"
    out.mkdir(exist_ok=True)
    mesh = roof_mesh(report)
    write_obj(report, mesh, out / "model.obj")
    write_ply(report, mesh, out / "model.ply", out / "edges.ply")
    preview(report, mesh, xyz, out / "model_preview.png")
    print(
        f"{args.house_id}: {len(mesh['faces'])} triangles over {len(set(mesh['facet_of_face'].tolist()))} facets, "
        f"mesh area {mesh_area(mesh):.1f} m² (report {report['totals']['roof_area_m2']} m²) -> {out}"
    )


if __name__ == "__main__":
    main()
