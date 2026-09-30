import numpy as np
import pytest
from shapely.geometry import Point, box

from conftest import DENSITY, OVERHANG
from roofrecon.isolate import (
    buffer_sweep,
    classification_histogram,
    crop_mask,
    density_stats,
    recenter,
)


def building_xy(scene):
    return scene["xyz"][scene["classification"] == 6, :2]


def test_classification_histogram(scene):
    h = classification_histogram(scene["classification"])
    assert set(h) == {2, 5, 6}
    assert sum(h.values()) == len(scene["xyz"])


def test_wall_footprint_clips_the_eaves(scene):
    """Cropping at the walls loses the overhang: the bias stage 5 has to account for."""
    n = crop_mask(building_xy(scene), scene["footprint"], 0.0).sum()
    assert n < scene["n_house"]
    walls_area = scene["footprint"].area
    roof_plan_area = scene["footprint"].buffer(OVERHANG, join_style="mitre").area
    assert n / scene["n_house"] == pytest.approx(walls_area / roof_plan_area, rel=0.05)


def test_buffer_past_eave_gets_whole_house_only(scene):
    n = crop_mask(building_xy(scene), scene["footprint"], OVERHANG + 0.1).sum()
    assert n == scene["n_house"]


def test_sweep_plateaus_then_rises_at_neighbour(scene):
    rows = buffer_sweep(building_xy(scene), scene["footprint"], [0.0, 0.5, 1.0, 2.0, 3.0])
    by_b = {r["buffer_m"]: r for r in rows}
    assert by_b[1.0]["added"] == 0 and by_b[2.0]["added"] == 0  # plateau: nothing between eave and neighbour
    # Neighbour's eave starts at NEIGHBOUR_GAP - OVERHANG = 2.6 m from our wall.
    assert by_b[3.0]["added"] > 0


def test_mitre_buffer_keeps_square_corners():
    fp = box(0, 0, 10, 8)
    corner_eave = np.array([[10.4, 8.4]])  # diagonal distance 0.566 m from the wall corner
    assert crop_mask(corner_eave, fp, 0.5)[0]
    assert not fp.buffer(0.5).contains(Point(10.4, 8.4))  # round join would drop it


def test_recenter_rescues_float32_precision(scene):
    xyz = scene["xyz"]
    local, offset = recenter(xyz)
    np.testing.assert_allclose(local + offset, xyz)
    assert np.all(offset == np.floor(offset))

    err_raw = np.abs(xyz[:, :2].astype(np.float32).astype(np.float64) - xyz[:, :2]).max()
    err_local = np.abs(local[:, :2].astype(np.float32).astype(np.float64) - local[:, :2]).max()
    assert err_raw > 5e-3  # ~cm-level rounding at y ~ 446 000 m
    assert err_local < 1e-5


def test_local_density_recovers_true_density(scene):
    xy = building_xy(scene)
    xy = xy[crop_mask(xy, scene["footprint"], 0.0)]
    d = density_stats(xy, scene["footprint"].area)
    assert d["nominal_pts_per_m2"] == pytest.approx(DENSITY, rel=0.05)
    # Edge points see half a disc, so the median sits slightly low; 15% on a 10 x 8 m roof.
    assert d["local_pts_per_m2"] == pytest.approx(DENSITY, rel=0.15)


def test_nn_spacing_on_a_grid():
    g = np.mgrid[0:10:0.25, 0:10:0.25].reshape(2, -1).T
    d = density_stats(g, 100.0)
    assert d["median_nn_spacing_m"] == pytest.approx(0.25)
    assert d["duplicate_xy_fraction"] == 0.0
