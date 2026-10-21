"""The 3D roof model: report facets as a triangle mesh, report edges as coloured lines.

Built from report.json, so the model and the measurement report always describe the same
roof. Coordinates stay local (stage 1 recentering): many viewers store vertices as float32,
which at Dutch map coordinates (y ~ 450 000 m) only resolves ~3 cm.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon

# Validated categorical palette (reference, light mode), in fixed order: one colour per facet,
# cycling only past eight facets. Edge colours follow stage 3's edge-type colours.
FACET_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
EDGE_COLORS = {"ridge": "#2a78d6", "hip": "#eb6834", "valley": "#1baf7a", "eave": "#eda100", "rake": "#e87ba4",
               "break": "#008300", "top": "#6250d6", "perimeter": "#e34948", "unclassified": "#52514e"}


def _rgb(hex_color: str) -> np.ndarray:
    return np.array([int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def _plane_frame(v: np.ndarray):
    """Origin and two in-plane axes for a planar ring of 3D points."""
    c = v.mean(axis=0)
    _, _, vt = np.linalg.svd(v - c)
    return c, vt[0], vt[1]


def triangulate_facet(vertices: np.ndarray, holes: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Triangles covering a planar facet polygon minus its holes: (points (K, 3), faces (T, 3)).

    Works in the facet's own plane: project to 2D, constrained Delaunay triangulation (which
    keeps the polygon's edges and holes), lift back. Only triangles inside the polygon are kept.
    """
    c, a, b = _plane_frame(vertices)
    to2 = lambda p: np.column_stack([(p - c) @ a, (p - c) @ b])  # noqa: E731
    shape = shapely.make_valid(Polygon(to2(vertices)))
    if holes:
        shape = shape.difference(shapely.union_all([shapely.make_valid(Polygon(to2(h))) for h in holes]))
    if shape.is_empty or shape.area <= 0:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=int)
    tris = shapely.get_parts(shapely.constrained_delaunay_triangles(shape))
    tris = [t for t in tris if shape.buffer(1e-6).contains(t.centroid)]
    pts, faces, index = [], [], {}
    for t in tris:
        face = []
        for xy in np.asarray(t.exterior.coords)[:3]:
            key = (round(xy[0], 6), round(xy[1], 6))
            if key not in index:
                index[key] = len(pts)
                pts.append(c + xy[0] * a + xy[1] * b)
            face.append(index[key])
        faces.append(face)
    return np.array(pts), np.array(faces, dtype=int)


def roof_mesh(report: dict) -> dict:
    """All counted facets as one mesh, with a facet id and colour per triangle."""
    verts, faces, facet_of_face = [], [], []
    for f in report["facets"]:
        if f["area_m2"] <= 0:
            continue
        p, t = triangulate_facet(np.asarray(f["vertices_local_m"]), [np.asarray(h) for h in f.get("holes_local_m", [])])
        if len(t) == 0:
            continue
        faces.append(t + sum(len(v) for v in verts))
        verts.append(p)
        facet_of_face += [f["id"]] * len(t)
    v = np.vstack(verts) if verts else np.zeros((0, 3))
    fc = np.vstack(faces) if faces else np.zeros((0, 3), dtype=int)
    return {"vertices": v, "faces": fc, "facet_of_face": np.array(facet_of_face, dtype=int)}


def mesh_area(mesh: dict) -> float:
    v, f = mesh["vertices"], mesh["faces"]
    if len(f) == 0:
        return 0.0
    return float(0.5 * np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1).sum())


def write_obj(report: dict, mesh: dict, path: Path) -> None:
    """OBJ + MTL: one group and material per facet, so viewers can pick and colour facets."""
    path = Path(path)
    mtl = path.with_suffix(".mtl")
    ids = [f["id"] for f in report["facets"] if f["id"] in set(mesh["facet_of_face"].tolist())]
    with mtl.open("w") as fh:
        for i, fid in enumerate(ids):
            r, g, b = _rgb(FACET_COLORS[i % len(FACET_COLORS)])
            fh.write(f"newmtl facet_{fid}\nKd {r:.4f} {g:.4f} {b:.4f}\nKa 0 0 0\nd 1.0\n\n")
    off = report.get("source", {}).get("offset_xyz")
    with path.open("w") as fh:
        fh.write(f"# roofrecon roof model: {report.get('house_id')}\n")
        fh.write(f"# local coordinates in metres; add offset {off} for {report.get('source', {}).get('crs')}\n")
        fh.write(f"mtllib {mtl.name}\n")
        for x, y, z in mesh["vertices"]:
            fh.write(f"v {x:.4f} {y:.4f} {z:.4f}\n")
        facets = {f["id"]: f for f in report["facets"]}
        for fid in ids:
            f = facets[fid]
            fh.write(f"g facet_{fid}_{f['facing'] or 'flat'}_{f['pitch_deg']}deg\nusemtl facet_{fid}\n")
            for a, b, c in mesh["faces"][mesh["facet_of_face"] == fid]:
                fh.write(f"f {a + 1} {b + 1} {c + 1}\n")


def write_ply(report: dict, mesh: dict, mesh_path: Path, edges_path: Path) -> None:
    """Mesh with per-vertex facet colours, and the edges as a coloured line set (Open3D)."""
    import open3d as o3d

    ids = [f["id"] for f in report["facets"] if f["id"] in set(mesh["facet_of_face"].tolist())]
    colour_of = {fid: _rgb(FACET_COLORS[i % len(FACET_COLORS)]) for i, fid in enumerate(ids)}
    # Duplicate vertices per face so each facet keeps its own colour at shared corners.
    v = mesh["vertices"][mesh["faces"].reshape(-1)]
    f = np.arange(len(v)).reshape(-1, 3)
    cols = np.repeat([colour_of[fid] for fid in mesh["facet_of_face"]], 3, axis=0)
    m = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v), o3d.utility.Vector3iVector(f))
    m.vertex_colors = o3d.utility.Vector3dVector(cols)
    m.compute_vertex_normals()
    o3d.io.write_triangle_mesh(str(mesh_path), m, write_ascii=False)

    pts, lines, lcols = [], [], []
    for e in report["edges"]:
        lines.append([len(pts), len(pts) + 1])
        pts += [e["p0_local_m"], e["p1_local_m"]]
        lcols.append(_rgb(EDGE_COLORS.get(e["type"], "#52514e")))
    ls = o3d.geometry.LineSet(o3d.utility.Vector3dVector(np.array(pts, float)), o3d.utility.Vector2iVector(np.array(lines)))
    ls.colors = o3d.utility.Vector3dVector(np.array(lcols))
    o3d.io.write_line_set(str(edges_path), ls)
