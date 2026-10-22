"""Interactive view of a house's roof model over its laser points (Open3D).

    uv run python scripts/view_model.py --house-id delft_03 [--no-points]

Facets coloured, edges coloured by type (see stage 3's legend), laser points in grey.
"""

import argparse
from pathlib import Path

import open3d as o3d

ap = argparse.ArgumentParser()
ap.add_argument("--house-id", required=True)
ap.add_argument("--no-points", action="store_true")
ap.add_argument("--root", type=Path, default=Path("data/interim"))
args = ap.parse_args()

house = args.root / args.house_id
geoms = [
    o3d.io.read_triangle_mesh(str(house / "model" / "model.ply")),
    o3d.io.read_line_set(str(house / "model" / "edges.ply")),
]
if not args.no_points:
    pcd = o3d.io.read_point_cloud(str(house / "roof.ply"))
    pcd.paint_uniform_color([0.32, 0.32, 0.31])
    geoms.append(pcd)
o3d.visualization.draw_geometries(geoms, window_name=f"{args.house_id} roof model", mesh_show_back_face=True)
