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


def region_grow(
    xyz: np.ndarray,
    normals: np.ndarray,
    variation: np.ndarray,
    grow_radius: float = 0.3,
    angle_deg: float = 8.0,
    threshold: float = 0.08,
    min_points: int = 30,
    seed_variation_factor: float = 2.0,
) -> Segmentation:
    """Flood-fill facets through neighbouring points whose normals agree.

    Seeds are taken flattest-first (lowest surface variation), so regions start in facet
    interiors, not on ridges. A neighbour joins a region if
      - its normal is within angle_deg of the region's plane normal, and
      - it lies within `threshold` of the region's plane.
    The region's plane is refitted each time the region doubles in size, so it starts
    from the seed's noisy normal but settles on the facet's.

    Comparing against the region's plane, not the neighbour that reached the point,
    stops slow drift around curved or gently bending surfaces. The normal test is what
    RANSAC lacks: points on a 15 degree slope have normals 15 degrees off vertical, so a
    horizontal region can't swallow them however close they are in height.

    Only points with surface variation <= seed_variation_factor x the median may seed.
    Along ridges and hips the normals are a blend of both facets, and that blend is
    consistent along the whole line, so a seed there grows a thin fake plane with an
    in-between pitch and facing. Band points can still join a real facet's region.

    Regions under min_points are released: their points can still join a later region,
    but never seed one.
    """
    n = len(xyz)
    tree = cKDTree(xyz)
    neighbours = tree.query_ball_point(xyz, grow_radius)
    cos_max = np.cos(np.radians(angle_deg))
    valid = ~np.isnan(normals[:, 0])
    order = np.argsort(np.where(valid, variation, np.inf))
    seed_ok = valid & (variation <= seed_variation_factor * np.nanmedian(variation))

    labels = np.full(n, -1)
    tried = np.zeros(n, dtype=bool)
    planes: list[Plane] = []
    for seed in order:
        if not seed_ok[seed] or labels[seed] != -1 or tried[seed]:
            continue
        k = len(planes)
        region = [seed]
        labels[seed] = k
        normal, origin = normals[seed], xyz[seed]
        next_refit = 10
        queue = [seed]
        while queue:
            p = queue.pop()
            for q in neighbours[p]:
                if labels[q] != -1 or not valid[q]:
                    continue
                if abs(normals[q] @ normal) < cos_max or abs((xyz[q] - origin) @ normal) >= threshold:
                    continue
                labels[q] = k
                region.append(q)
                queue.append(q)
            if len(region) >= next_refit:
                fitted = fit_plane(xyz[region])
                normal, origin = fitted.normal, fitted.centroid
                next_refit *= 2
        tried[region] = True
        if len(region) < min_points:
            labels[region] = -1
            continue
        planes.append(fit_plane(xyz[region]))
    return Segmentation(labels, planes)


def absorb_unassigned(xyz: np.ndarray, seg: Segmentation, threshold: float, radius: float = 0.3, rounds: int = 10) -> Segmentation:
    """Give unassigned points to an adjacent plane they lie on.

    Region growing leaves the band along ridges unassigned, because normals there blend
    both facets. A point joins a plane if a neighbour within `radius` already belongs to
    it and the point is within `threshold` of it; if several qualify, the nearest plane
    wins. Repeating lets the facets grow into the band from both sides.
    """
    labels = seg.labels.copy()
    tree = cKDTree(xyz)
    for _ in range(rounds):
        todo = np.flatnonzero(labels == -1)
        if len(todo) == 0:
            break
        changed = False
        for i, nbrs in zip(todo, tree.query_ball_point(xyz[todo], radius)):
            cands = {labels[j] for j in nbrs if labels[j] >= 0}
            if not cands:
                continue
            dists = {k: abs(seg.planes[k].distance(xyz[i:i + 1])[0]) for k in cands}
            best = min(dists, key=dists.get)
            if dists[best] < threshold:
                labels[i] = best
                changed = True
        if not changed:
            break
    planes = [fit_plane(xyz[labels == k]) for k in range(len(seg.planes))]
    return Segmentation(labels, planes)


def segment_region_growing(
    xyz: np.ndarray,
    normal_radius: float = 0.5,
    angle_deg: float = 8.0,
    threshold: float = 0.08,
    grow_radius: float = 0.3,
    min_points: int = 30,
) -> Segmentation:
    """Normals -> region growing -> merge -> absorb the ridge bands -> reassign -> merge -> reassign.

    The second merge catches pieces of one facet that were separated by the unassigned
    band. They only touch, and only fit one plane cleanly, once the band is filled and
    ridge points have gone to their nearest plane.
    """
    from roofrecon.normals import estimate_normals

    normals, variation, _ = estimate_normals(xyz, normal_radius)
    seg = region_grow(xyz, normals, variation, grow_radius, angle_deg, threshold, min_points)
    seg = merge_coplanar(xyz, seg, touch_radius=grow_radius)
    seg = absorb_unassigned(xyz, seg, threshold, grow_radius)
    seg = reassign_to_nearest_plane(xyz, seg, threshold)
    seg = merge_coplanar(xyz, seg, touch_radius=grow_radius)
    return reassign_to_nearest_plane(xyz, seg, threshold)


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
