"""Open a .ply in the Open3D viewer:  uv run python scripts/view.py data/interim/<house-id>/roof.ply"""

import sys

import open3d as o3d

pcd = o3d.io.read_point_cloud(sys.argv[1])
print(pcd)
axes = o3d.geometry.TriangleMesh.create_coordinate_frame(size=2.0)  # 2 m axes: red x, green y, blue z
o3d.visualization.draw_geometries([pcd, axes], window_name=sys.argv[1])
