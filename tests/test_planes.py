import numpy as np
import pytest

from roofrecon.planes import (
    azimuth_deg,
    fit_plane,
    hull_area,
    pitch_deg,
    plan_hull_area,
    plane_axes,
)
from synthetic_roofs import gable, hip


def test_fit_recovers_exact_plane_without_noise():
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 10, (200, 2))
    z = 2.0 + 0.5 * xy[:, 0] - 0.25 * xy[:, 1]  # z = 2 + 0.5x - 0.25y
    p = fit_plane(np.column_stack([xy, z]))
    expected = np.array([-0.5, 0.25, 1.0]) / np.linalg.norm([-0.5, 0.25, 1.0])
    np.testing.assert_allclose(p.normal, expected, atol=1e-9)
    assert p.rms < 1e-9


def test_normal_always_points_up():
    xyz = np.array([[0, 0, 0], [1, 0, 1], [0, 1, 0], [1, 1, 1]], dtype=float)
    assert fit_plane(xyz).normal[2] > 0
    assert fit_plane(xyz[::-1]).normal[2] > 0


def test_flat_plane_has_zero_pitch_and_no_azimuth():
    xyz = np.array([[0, 0, 3], [4, 0, 3], [0, 4, 3], [4, 4, 3]], dtype=float)
    p = fit_plane(xyz)
    assert pitch_deg(p.normal) == pytest.approx(0.0, abs=1e-9)
    assert azimuth_deg(p.normal) is None


@pytest.mark.parametrize("pitch", [15.0, 35.0, 50.0])
def test_gable_pitch_recovered_within_half_degree(pitch):
    roof = gable(pitch_deg=pitch)
    for k in (0, 1):
        p = fit_plane(roof.xyz[roof.labels == k])
        assert pitch_deg(p.normal) == pytest.approx(pitch, abs=0.5)
        # The 3 cm noise is vertical; the fit measures it perpendicular to the plane, which
        # shrinks it by cos(pitch). Same for real LiDAR, whose error is mostly along the beam.
        assert p.rms == pytest.approx(0.03 * np.cos(np.radians(pitch)), rel=0.15)


def test_azimuth_is_downhill_direction():
    roof = gable()
    south = fit_plane(roof.xyz[roof.labels == 0])  # slopes down towards y = 0
    north = fit_plane(roof.xyz[roof.labels == 1])
    assert azimuth_deg(south.normal) == pytest.approx(180.0, abs=1.0)
    north_az = azimuth_deg(north.normal)
    assert min(north_az, 360 - north_az) < 1.0  # within 1 degree of north, either side of 0/360


def test_plane_axes_are_orthonormal_and_in_plane():
    p = fit_plane(gable().xyz[gable().labels == 0])
    u, v = plane_axes(p.normal)
    for a in (u, v):
        assert np.linalg.norm(a) == pytest.approx(1.0)
        assert a @ p.normal == pytest.approx(0.0, abs=1e-12)
    assert u @ v == pytest.approx(0.0, abs=1e-12)
    assert u[2] == pytest.approx(0.0, abs=1e-12)  # u is horizontal
    assert v[2] > 0  # v points up the slope


@pytest.mark.parametrize("make", [gable, hip])
def test_true_area_matches_ground_truth(make):
    roof = make()
    for k, truth in enumerate(roof.true_areas):
        pts = roof.xyz[roof.labels == k]
        area = hull_area(pts, fit_plane(pts))
        # The hull sits slightly inside the real edges, so it reads a little low, never high
        # (all synthetic facets are convex).
        assert area == pytest.approx(truth, rel=0.03)
        assert area < truth


def test_true_area_is_plan_area_over_cos_pitch():
    roof = gable(pitch_deg=40.0)
    pts = roof.xyz[roof.labels == 0]
    p = fit_plane(pts)
    ratio = hull_area(pts, p) / plan_hull_area(pts)
    assert ratio == pytest.approx(1 / np.cos(np.radians(40.0)), rel=0.01)
