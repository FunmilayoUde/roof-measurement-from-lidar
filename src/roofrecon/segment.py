"""Splitting roof points into planar facets."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from roofrecon.planes import Plane, fit_plane


@dataclass
class Segmentation:
    labels: np.ndarray  # (N,) facet index per point, -1 = unassigned
    planes: list[Plane] = field(default_factory=list)

    @property
    def coverage(self) -> float:
        """Fraction of points assigned to some plane."""
        return float(np.mean(self.labels >= 0)) if len(self.labels) else 0.0


def ransac_iterations(inlier_fraction: float, confidence: float = 0.999, sample_size: int = 3) -> int:
    """How many random samples are needed to draw at least one all-inlier sample.

    One sample of 3 points is all-inlier with probability w^3. The chance that N samples
    all miss is (1 - w^3)^N; setting that to 1 - confidence gives
    N = log(1 - confidence) / log(1 - w^3). With w = 0.25 (a facet holding a quarter of
    the remaining points) and 99.9% confidence, N is about 440.
    """
    p_good = inlier_fraction**sample_size
    if p_good >= 1.0:
        return 1
    return int(np.ceil(np.log(1 - confidence) / np.log(1 - p_good)))


def ransac_plane(
    xyz: np.ndarray, threshold: float, iterations: int, rng: np.random.Generator
) -> tuple[Plane | None, np.ndarray]:
    """Best plane by inlier count, refined by least squares on its inliers.

    The 3-point plane from the winning sample is exact through three noisy points, so it
    is tilted by their noise. Refitting on all inliers averages the noise out.
    """
    n = len(xyz)
    best_count, best_mask = 0, None
    for _ in range(iterations):
        a, b, c = xyz[rng.choice(n, 3, replace=False)]
        normal = np.cross(b - a, c - a)
        norm = np.linalg.norm(normal)
        if norm < 1e-9:  # three nearly collinear points don't define a plane
            continue
        normal /= norm
        mask = np.abs((xyz - a) @ normal) < threshold
        count = int(mask.sum())
        if count > best_count:
            best_count, best_mask = count, mask
    if best_mask is None or best_count < 3:
        return None, np.zeros(n, dtype=bool)
    plane = fit_plane(xyz[best_mask])
    return plane, np.abs(plane.distance(xyz)) < threshold


def largest_component(xyz: np.ndarray, mask: np.ndarray, radius: float) -> np.ndarray:
    """Keep only the largest spatially connected group of the masked points.

    Two points are connected if they are within `radius` (in 3D). RANSAC judges points
    only by distance to an infinite plane, so it happily collects separate patches that
    happen to be coplanar: two flat roofs at the same height, or a strip of the opposite
    facet near a ridge. Real facets are contiguous, so the extra patches go back in the
    pool for later planes.
    """
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return mask
    pairs = cKDTree(xyz[idx]).query_pairs(radius, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(idx), len(idx)))
    _, comp = connected_components(graph, directed=False)
    keep = idx[comp == np.bincount(comp).argmax()]
    out = np.zeros_like(mask)
    out[keep] = True
    return out


def sequential_ransac(
    xyz: np.ndarray,
    threshold: float = 0.08,
    min_points: int = 30,
    connect_radius: float | None = 0.3,
    iterations: int = 500,
    max_planes: int = 30,
    seed: int = 0,
) -> Segmentation:
    """Find planes one at a time: fit the best plane to the remaining points, remove its
    inliers, repeat until the best plane has fewer than min_points.

    - threshold: max point-to-plane distance (perpendicular, metres) to count as inlier.
    - min_points: smallest plane worth keeping. Think of it as min area x density.
    - connect_radius: if set, keep only the largest connected patch of each plane's
      inliers. Must exceed the typical point spacing or real facets fall apart.
    """
    rng = np.random.default_rng(seed)
    labels = np.full(len(xyz), -1)
    planes: list[Plane] = []
    remaining = np.arange(len(xyz))
    while len(remaining) >= max(min_points, 3) and len(planes) < max_planes:
        sub = xyz[remaining]
        plane, mask = ransac_plane(sub, threshold, iterations, rng)
        if plane is None:
            break
        if connect_radius is not None:
            mask = largest_component(sub, mask, connect_radius)
        if mask.sum() < min_points:
            break
        plane = fit_plane(sub[mask])
        labels[remaining[mask]] = len(planes)
        planes.append(plane)
        remaining = remaining[~mask]
    return Segmentation(labels, planes)


def segment_ransac(
    xyz: np.ndarray,
    threshold: float = 0.08,
    min_points: int = 30,
    connect_radius: float | None = 0.3,
    iterations: int = 500,
    seed: int = 0,
) -> Segmentation:
    """Sequential RANSAC, then clean-up: reassign -> merge -> reassign.

    The order matters. Earlier planes steal strips of later facets along shared edges;
    at a narrow hip-end tip two stolen strips can meet and cut the tip off, which the
    connectivity check then turns into its own plane. Reassigning first gives the strips
    back, which reconnects the tip so the merge can see it touches its facet.
    """
    seg = sequential_ransac(xyz, threshold, min_points, connect_radius, iterations, seed=seed)
    seg = reassign_to_nearest_plane(xyz, seg, threshold)
    seg = merge_coplanar(xyz, seg, touch_radius=connect_radius or 0.3)
    return reassign_to_nearest_plane(xyz, seg, threshold)


def merge_coplanar(
    xyz: np.ndarray, seg: Segmentation, max_angle_deg: float = 5.0, rms_factor: float = 1.5, touch_radius: float = 0.3
) -> Segmentation:
    """Merge pairs of planes that are really one facet.

    Greedy RANSAC leaves scraps: once a facet's plane takes most of its points, the
    leftovers (often a strip along a ridge or hip) can form a second, almost identical
    plane. Two planes are merged when their normals differ by less than max_angle_deg,
    their point sets touch (within touch_radius), and one plane fits the union nearly as
    tightly as the two did separately (rms <= rms_factor x the worse of the two).
    The rms test is what keeps two parallel but offset planes, like a roof and a flat
    dormer top, apart.
    """
    labels = seg.labels.copy()
    planes = list(seg.planes)
    cos_max = np.cos(np.radians(max_angle_deg))
    merged = True
    while merged:
        merged = False
        trees = {k: cKDTree(xyz[labels == k]) for k in range(len(planes)) if planes[k] is not None}
        for a in list(trees):
            for b in list(trees):
                if b <= a or planes[a] is None or planes[b] is None:
                    continue
                if abs(planes[a].normal @ planes[b].normal) < cos_max:
                    continue
                if trees[a].count_neighbors(trees[b], touch_radius) == 0:
                    continue
                union = fit_plane(xyz[(labels == a) | (labels == b)])
                if union.rms > rms_factor * max(planes[a].rms, planes[b].rms):
                    continue
                labels[labels == b] = a
                planes[a], planes[b] = union, None
                merged = True
                break
            if merged:
                break
    # Renumber so labels stay 0..K-1.
    keep = [k for k, p in enumerate(planes) if p is not None]
    remap = np.full(len(planes), -1)
    remap[keep] = np.arange(len(keep))
    new_labels = np.where(labels >= 0, remap[np.maximum(labels, 0)], -1)
    return Segmentation(new_labels, [planes[k] for k in keep])


def reassign_to_nearest_plane(xyz: np.ndarray, seg: Segmentation, threshold: float) -> Segmentation:
    """Give each point to the closest plane it is an inlier of.

    Sequential RANSAC is greedy: the first plane keeps every point within `threshold`,
    including a strip of the next facet along their shared ridge. Both planes pass
    through the ridge line, so near it a point is closer to its own facet's plane.
    The strip is about threshold / sin(angle between facets) wide, so it is widest
    on shallow roofs. Unassigned points stay unassigned.
    """
    if not seg.planes:
        return seg
    dist = np.abs(np.column_stack([p.distance(xyz) for p in seg.planes]))
    nearest = dist.argmin(axis=1)
    labels = np.where((seg.labels >= 0) & (dist.min(axis=1) < threshold), nearest, seg.labels)
    planes = [fit_plane(xyz[labels == k]) if (labels == k).sum() >= 3 else p for k, p in enumerate(seg.planes)]
    return Segmentation(labels, planes)
