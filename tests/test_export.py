import numpy as np
import pytest

from roofrecon.export import mesh_area, roof_mesh, triangulate_facet, write_obj
from roofrecon.measure import measurement_report
from roofrecon.pipeline import build_roof
from synthetic_roofs import dormer, gable, hip


def test_triangulated_square_with_hole():
    sq = np.array([[0, 0, 1], [4, 0, 1], [4, 4, 1], [0, 4, 1]], float)
    hole = np.array([[1, 1, 1], [2, 1, 1], [2, 2, 1], [1, 2, 1]], float)
    p, t = triangulate_facet(sq, [hole])
    area = 0.5 * np.linalg.norm(np.cross(p[t[:, 1]] - p[t[:, 0]], p[t[:, 2]] - p[t[:, 0]]), axis=1).sum()
    assert area == pytest.approx(15.0)
    assert np.allclose(p[:, 2], 1.0)  # stays in the facet's plane


def test_sloped_facet_area_is_true_area():
    # A 4 x 4 square tilted 40 degrees about x: true area 16, plan area 16 cos 40.
    v = np.array(
        [
            [0, 0, 0],
            [4, 0, 0],
            [4, 4 * np.cos(np.radians(40)), 4 * np.sin(np.radians(40))],
            [0, 4 * np.cos(np.radians(40)), 4 * np.sin(np.radians(40))],
        ],
        float,
    )
    p, f = triangulate_facet(v, [])
    area = 0.5 * np.linalg.norm(np.cross(p[f[:, 1]] - p[f[:, 0]], p[f[:, 2]] - p[f[:, 0]]), axis=1).sum()
    assert area == pytest.approx(16.0)


@pytest.mark.parametrize("make", [gable, hip, dormer])
def test_mesh_area_matches_report(make, tmp_path):
    rep = measurement_report(build_roof(make().xyz, 33.0), {"house_id": make.__name__})
    mesh = roof_mesh(rep)
    # Traced facets trimmed in the report keep their full outline in the mesh, so compare
    # with the untrimmed sum.
    untrimmed = sum(f["area_m2"] + (f["area_not_counted_m2"] or 0) for f in rep["facets"] if f["area_m2"] > 0)
    assert mesh_area(mesh) == pytest.approx(untrimmed, rel=1e-3)
    write_obj(rep, mesh, tmp_path / "m.obj")
    text = (tmp_path / "m.obj").read_text()
    assert text.count("\nf ") == len(mesh["faces"]) and "mtllib m.mtl" in text
