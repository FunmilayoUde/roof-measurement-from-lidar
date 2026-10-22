import numpy as np

from roofrecon.normals import estimate_normals
from synthetic_roofs import gable


def true_normals(roof):
    t = np.radians(roof.pitch_deg)
    return np.where(roof.labels[:, None] == 0, [0, -np.sin(t), np.cos(t)], [0, np.sin(t), np.cos(t)])


def angle_error_deg(n, truth):
    return np.degrees(np.arccos(np.clip(np.abs(np.sum(n * truth, axis=1)), 0, 1)))


def test_interior_normals_accurate_at_half_metre():
    roof = gable()
    n, _, _ = estimate_normals(roof.xyz, 0.5)
    far = np.abs(roof.xyz[:, 1] - 4.0) > 1.0  # more than 1 m from the ridge
    assert np.median(angle_error_deg(n[far], true_normals(roof)[far])) < 2.0


def test_larger_radius_less_noise_but_wider_ridge_blur():
    roof = gable()
    d = np.abs(roof.xyz[:, 1] - 4.0)
    errs = {}
    for r in (0.3, 1.0):
        n, _, _ = estimate_normals(roof.xyz, r, max_nn=300)
        e = angle_error_deg(n, true_normals(roof))
        # nanmedian: at r = 0.3 some points have < 3 neighbours and a NaN normal
        errs[r] = (np.nanmedian(e[d > 1.5]), np.nanmedian(e[(d > 0.3) & (d < 0.4)]))
    assert errs[1.0][0] < errs[0.3][0]  # far from ridge: big radius wins
    assert errs[1.0][1] > errs[0.3][1]  # 0.3-0.4 m from ridge: big radius blurs


def test_normals_point_up_and_variation_peaks_at_ridge():
    roof = gable()
    n, var, _ = estimate_normals(roof.xyz, 0.5)
    assert np.all(n[:, 2] > 0)
    d = np.abs(roof.xyz[:, 1] - 4.0)
    assert np.median(var[d < 0.1]) > 3 * np.median(var[d > 1.0])


def test_too_few_neighbours_gives_nan():
    xyz = np.array([[0, 0, 0], [5, 0, 0], [0, 5, 0], [5, 5, 0.0]])
    n, var, cnt = estimate_normals(xyz, 1.0)
    assert np.all(cnt == 1) and np.all(np.isnan(n)) and np.all(np.isnan(var))
