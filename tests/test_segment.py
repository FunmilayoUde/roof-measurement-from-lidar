import numpy as np
import pytest

from roofrecon.planes import pitch_deg
from roofrecon.normals import estimate_normals
from roofrecon.segment import (
    merge_coplanar,
    ransac_iterations,
    region_grow,
    segment_ransac,
    segment_region_growing,
    sequential_ransac,
)
from synthetic_roofs import gable, hip


def mislabelled_fraction(truth: np.ndarray, labels: np.ndarray) -> float:
    """Points not in their segment's majority facet, over all points."""
    wrong = 0
    for k in np.unique(labels[labels >= 0]):
        m = labels == k
        wrong += m.sum() - np.bincount(truth[m]).max()
    return wrong / len(truth)


def slab(rng, x0, z, n_per_m2=33, length=5.0, width=4.0, noise=0.03):
    n = rng.poisson(n_per_m2 * length * width)
    return np.column_stack([rng.uniform(x0, x0 + length, n), rng.uniform(0, width, n), z + rng.normal(0, noise, n)])


def test_ransac_iterations_formula():
    assert ransac_iterations(0.25) == 439
    assert ransac_iterations(0.5) == 52
    assert ransac_iterations(1.0) == 1


def test_connectivity_splits_separate_coplanar_patches():
    rng = np.random.default_rng(0)
    xyz = np.vstack([slab(rng, 0, 3.0), slab(rng, 8, 3.0)])  # same height, 3 m apart
    assert len(sequential_ransac(xyz, connect_radius=None).planes) == 1
    assert len(sequential_ransac(xyz, connect_radius=0.3).planes) == 2


def test_merge_keeps_parallel_offset_planes_apart():
    """A 20 cm step between two flat roofs: same normal, touching, but not one plane."""
    rng = np.random.default_rng(0)
    xyz = np.vstack([slab(rng, 0, 3.0), slab(rng, 5, 3.2)])
    seg = sequential_ransac(xyz, threshold=0.08)
    assert len(merge_coplanar(xyz, seg).planes) == 2


def test_gable_two_planes_correct_pitch():
    roof = gable()
    seg = segment_ransac(roof.xyz, threshold=0.10)
    assert len(seg.planes) == 2
    assert seg.coverage > 0.99
    for p in seg.planes:
        assert pitch_deg(p.normal) == pytest.approx(35.0, abs=0.5)
    assert mislabelled_fraction(roof.labels, seg.labels) < 0.01


def test_cleanup_undoes_ridge_stealing_on_shallow_gable():
    """At 10 degrees and a 20 cm threshold, the first plane takes a wide strip of the
    other facet along the ridge; reassignment gives it back."""
    roof = gable(pitch_deg=10.0)
    raw = sequential_ransac(roof.xyz, threshold=0.20)
    clean = segment_ransac(roof.xyz, threshold=0.20)
    assert mislabelled_fraction(roof.labels, raw.labels) > 0.10
    assert mislabelled_fraction(roof.labels, clean.labels) < 0.04
    for p in clean.planes:
        assert pitch_deg(p.normal) == pytest.approx(10.0, abs=0.5)


@pytest.mark.parametrize("seed", range(5))
def test_hip_four_planes_after_cleanup(seed):
    roof = hip(seed=seed)
    seg = segment_ransac(roof.xyz, threshold=0.10, seed=seed)
    assert len(seg.planes) == 4
    assert mislabelled_fraction(roof.labels, seg.labels) < 0.025


def test_small_object_left_unassigned():
    """A chimney-sized cluster (20 points) is below min_points and stays unassigned."""
    roof = gable()
    rng = np.random.default_rng(1)
    chimney = np.column_stack([rng.uniform(4, 4.6, 20), rng.uniform(3.6, 4.2, 20), rng.uniform(8.0, 9.0, 20)])
    xyz = np.vstack([roof.xyz, chimney])
    seg = segment_ransac(xyz, threshold=0.10, min_points=30)
    assert np.all(seg.labels[-20:] == -1)
    assert len(seg.planes) == 2


def test_known_failure_horizontal_slab_on_shallow_hip():
    """Documents a RANSAC failure, not desired behaviour. On a 15 degree hip with a 20 cm
    threshold, a horizontal slab 40 cm thick through the ridge area holds more points than
    any real facet, so it wins. Distance alone can't reject it; per-point normals can.
    test_region_growing_avoids_horizontal_slab_on_shallow_hip covers the fix."""
    roof = hip(pitch_deg=15.0)
    seg = sequential_ransac(roof.xyz, threshold=0.20)
    assert pitch_deg(seg.planes[0].normal) < 5.0


def test_region_growing_avoids_horizontal_slab_on_shallow_hip():
    """The same case RANSAC fails above: normals tilted 15 degrees can't join a flat region."""
    roof = hip(pitch_deg=15.0)
    seg = segment_region_growing(roof.xyz, threshold=0.20)
    assert len(seg.planes) == 4
    for p in seg.planes:
        assert pitch_deg(p.normal) == pytest.approx(15.0, abs=0.5)
    assert mislabelled_fraction(roof.labels, seg.labels) < 0.05


def test_region_growing_gable():
    roof = gable()
    seg = segment_region_growing(roof.xyz, threshold=0.10)
    assert len(seg.planes) == 2
    assert seg.coverage > 0.99
    assert mislabelled_fraction(roof.labels, seg.labels) < 0.01


@pytest.mark.parametrize("seed", range(5))
def test_region_growing_hip(seed):
    roof = hip(seed=seed)
    seg = segment_region_growing(roof.xyz, threshold=0.10)
    assert len(seg.planes) == 4
    assert mislabelled_fraction(roof.labels, seg.labels) < 0.025


def test_ridge_band_left_by_growing_is_absorbed():
    roof = hip()
    n, var, _ = estimate_normals(roof.xyz, 0.5)
    grown = region_grow(roof.xyz, n, var, threshold=0.10)
    final = segment_region_growing(roof.xyz, threshold=0.10)
    assert grown.coverage < 0.90  # blended normals along ridges and hips don't join (88.5%)
    # 98.8%: absorption needs a labelled neighbour, so a few sparse points at the eave
    # corners whose only neighbours are also unassigned stay out.
    assert final.coverage > 0.98
