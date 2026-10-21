import numpy as np
import pytest

from roofrecon.planes import Plane, fit_plane
from roofrecon.segment import Segmentation, segment_ransac
from roofrecon.topology import internal_edges, plane_intersection, three_plane_point
from synthetic_roofs import cross_gable, gable, hip


def plane(normal, point) -> Plane:
    n = np.asarray(normal, float)
    return Plane(n / np.linalg.norm(n), np.asarray(point, float), 0.0)


def truth_segmentation(roof) -> Segmentation:
    k = roof.labels.max() + 1
    return Segmentation(roof.labels.copy(), [fit_plane(roof.xyz[roof.labels == i]) for i in range(k)])


def match(roof, edges, mapping):
    """Pair each found edge with the ground-truth edge between the same two facets."""
    out = []
    for e in edges:
        key = tuple(sorted(mapping[f] for f in e.facets))
        truth = [t for t in roof.edges if tuple(sorted(t["facets"])) == key]
        out.append((e, truth[0] if truth else None))
    return out


def endpoint_errors(e, t):
    """Distance of each found endpoint to the nearer true endpoint."""
    return [min(np.linalg.norm(p - t["p0"]), np.linalg.norm(p - t["p1"])) for p in (e.p0, e.p1)]


def test_plane_intersection_of_two_roof_slopes():
    # z = 5 + y and z = 13 - y meet along y = 4, z = 9, running along x.
    a = plane([0, -1, 1], [0, 0, 5])
    b = plane([0, 1, 1], [0, 8, 5])
    point, d = plane_intersection(a, b)
    assert abs(d[0]) == pytest.approx(1.0)
    assert point[1] == pytest.approx(4.0) and point[2] == pytest.approx(9.0)


def test_parallel_planes_have_no_intersection():
    assert plane_intersection(plane([0, 0, 1], [0, 0, 3]), plane([0, 0, 1], [0, 0, 3.2])) is None


def test_three_plane_point():
    q = three_plane_point(plane([1, 0, 0], [2, 0, 0]), plane([0, 1, 0], [0, 3, 0]), plane([0, 0, 1], [0, 0, 4]))
    np.testing.assert_allclose(q, [2, 3, 4])


@pytest.mark.parametrize("make", [gable, hip, cross_gable])
@pytest.mark.parametrize("use_ransac", [False, True])
def test_edge_types_and_count(make, use_ransac):
    roof = make()
    if use_ransac:
        seg = segment_ransac(roof.xyz, threshold=0.12, min_points=33)
        mapping = {k: np.bincount(roof.labels[seg.labels == k]).argmax() for k in range(len(seg.planes))}
    else:
        seg, mapping = truth_segmentation(roof), {k: k for k in range(roof.labels.max() + 1)}
    pairs = match(roof, internal_edges(roof.xyz, seg), mapping)
    assert len(pairs) == len(roof.edges)
    for e, t in pairs:
        assert t is not None, f"extra edge between {e.facets}"
        assert e.kind == t["type"]


@pytest.mark.parametrize("make", [hip, cross_gable])
def test_corner_snapped_endpoints_are_exact(make):
    """Where three facets meet, the endpoint is the three-plane intersection: exact to the
    plane-fit accuracy, not limited by point spacing."""
    roof = make()
    for e, t in match(roof, internal_edges(roof.xyz, truth_segmentation(roof)), {k: k for k in range(4)}):
        for err, snapped in zip(endpoint_errors(e, t), e.snapped):
            if snapped:
                assert err < 0.05


def test_hip_ridge_length_exact_when_both_ends_snap():
    roof = hip()
    edges = internal_edges(roof.xyz, truth_segmentation(roof))
    ridge = [e for e in edges if e.kind == "ridge"][0]
    assert ridge.snapped == (True, True)
    assert ridge.length == pytest.approx(4.0, abs=0.02)


@pytest.mark.parametrize("make", [gable, hip, cross_gable])
def test_unsnapped_ends_fall_short_never_long(make):
    """Ends without a three-plane corner come from where the contact points stop, so they
    can only be short. Step 2 replaces them with intersections with the outer edges."""
    roof = make()
    for e, t in match(roof, internal_edges(roof.xyz, truth_segmentation(roof)), {k: k for k in range(4)}):
        assert e.length <= t["length"] + 0.02


def test_mansard_curb_is_a_break():
    """A steep side meeting a flat top: only one facet slopes away, so neither ridge nor valley."""
    rng = np.random.default_rng(0)
    n = 3000
    x, y = rng.uniform(0, 10, n), rng.uniform(0, 4, n)
    t55 = np.tan(np.radians(55))
    z = np.where(y >= 2, 6.0, 6.0 - t55 * (2 - y)) + rng.normal(0, 0.02, n)
    labels = (y >= 2).astype(int)
    xyz = np.column_stack([x, y, z])
    seg = Segmentation(labels, [fit_plane(xyz[labels == k]) for k in (0, 1)])
    edges = internal_edges(xyz, seg)
    assert [e.kind for e in edges] == ["break"]


@pytest.mark.parametrize("seed", [2, 5, 9, 15])
def test_valley_beside_notched_facet(seed):
    """A's south slope wraps around the wing, so it has points on both sides of each valley
    line. Deciding its side from all its points (centroid) misread these seeds' valleys as
    breaks; the side must come from points near the edge."""
    roof = cross_gable(seed=seed)
    seg = segment_ransac(roof.xyz, threshold=0.12, min_points=33, seed=seed)
    mapping = {k: np.bincount(roof.labels[seg.labels == k]).argmax() for k in range(len(seg.planes))}
    for e, t in match(roof, internal_edges(roof.xyz, seg), mapping):
        assert e.kind == t["type"]


@pytest.mark.parametrize("seed", [8, 16])
def test_valley_short_contact_and_mixed_labels(seed):
    """Seed 8: contact stopped short, both ends were within snap distance of the same corner
    and the edge collapsed to zero length. Seed 16: labels right at the valley line mix, so
    the side test must look a little way back from the line."""
    roof = cross_gable(seed=seed)
    seg = segment_ransac(roof.xyz, threshold=0.12, min_points=33, seed=seed)
    mapping = {k: np.bincount(roof.labels[seg.labels == k]).argmax() for k in range(len(seg.planes))}
    for e, t in match(roof, internal_edges(roof.xyz, seg), mapping):
        assert e.kind == t["type"]
        assert e.length > 0.5
