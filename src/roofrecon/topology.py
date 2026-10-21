"""Roof topology: which facets meet, the lines where they meet, and what kind of edge each is."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from roofrecon.planes import Plane
from roofrecon.segment import Segmentation

HORIZONTAL_TOL_DEG = 10.0  # an edge within this of horizontal counts as horizontal (ridge, not hip)
FLAT_TOL = np.sin(np.radians(5.0))  # facet slope below ~5 degrees away from the edge counts as flat


@dataclass
class Edge:
    kind: str  # ridge | hip | valley | break
    facets: tuple[int, int]
    p0: np.ndarray
    p1: np.ndarray
    snapped: tuple[bool, bool] = (False, False)  # endpoint fixed by a three-plane corner?

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))


def plane_intersection(a: Plane, b: Plane, min_angle_deg: float = 2.0) -> tuple[np.ndarray, np.ndarray] | None:
    """The line where two planes meet, as (point, unit direction), or None if nearly parallel.

    The direction lies in both planes, so it is perpendicular to both normals: their cross
    product. For a point, solve the two plane equations plus a third that pins the point
    to the plane through a's centroid perpendicular to the line, which makes the system
    square and well conditioned.
    """
    d = np.cross(a.normal, b.normal)
    s = np.linalg.norm(d)
    if s < np.sin(np.radians(min_angle_deg)):
        return None
    d /= s
    A = np.vstack([a.normal, b.normal, d])
    rhs = np.array([a.normal @ a.centroid, b.normal @ b.centroid, d @ a.centroid])
    return np.linalg.solve(A, rhs), d


def three_plane_point(a: Plane, b: Plane, c: Plane, min_det: float = 0.05) -> np.ndarray | None:
    """The single point where three planes meet, or None if two are nearly parallel."""
    A = np.vstack([a.normal, b.normal, c.normal])
    if abs(np.linalg.det(A)) < min_det:
        return None
    return np.linalg.solve(A, np.array([a.normal @ a.centroid, b.normal @ b.centroid, c.normal @ c.centroid]))


def adjacency(xyz: np.ndarray, labels: np.ndarray, n_planes: int, radius: float = 0.3) -> dict[tuple[int, int], np.ndarray]:
    """Pairs of facets that touch, with their contact points.

    Contact points are the points of either facet that have a point of the other within
    `radius`. They sit along the shared edge, so they say where along the infinite
    intersection line the real edge is.
    """
    trees = {k: cKDTree(xyz[labels == k]) for k in range(n_planes) if (labels == k).any()}
    idx = {k: np.flatnonzero(labels == k) for k in trees}
    pairs = {}
    for a in trees:
        for b in trees:
            if b <= a:
                continue
            near_b = trees[b].query_ball_point(xyz[idx[a]], radius, return_length=True) > 0
            if not near_b.any():
                continue
            near_a = trees[a].query_ball_point(xyz[idx[b]], radius, return_length=True) > 0
            pairs[(a, b)] = np.concatenate([idx[a][near_b], idx[b][near_a]])
    return pairs


def slope_away(plane: Plane, point: np.ndarray, direction: np.ndarray, side_points: np.ndarray) -> float:
    """Vertical component of the unit vector that lies in `plane`, is perpendicular to the
    edge line, and points from the line into the facet (towards side_points).
    Negative: the facet falls away from the edge (convex side). Positive: it rises.

    side_points must be the facet's points *beside this edge*: alongside it, and a little
    way back from it. A notched facet (a main slope with a wing cutting into it) has points
    on both sides of the infinite edge line, so its overall centroid can be on the wrong
    side; and right at the line, labels from both facets mix. The median of the offsets
    decides the side, so a few mislabelled points can't flip it."""
    u = np.cross(plane.normal, direction)
    u /= np.linalg.norm(u)
    rel = side_points - point
    if np.median((rel - np.outer(rel @ direction, direction)) @ u) < 0:
        u = -u
    return float(u[2])


def beside_edge(pts: np.ndarray, point, direction, t0: float, t1: float, near: float = 0.2, far: float = 1.5) -> np.ndarray:
    """A facet's points alongside the edge segment and between `near` and `far` from the line.
    Falls back to the nearest ones if the strip is empty (a very narrow facet)."""
    rel = pts - point
    along = rel @ direction
    dist = np.linalg.norm(rel - np.outer(along, direction), axis=1)
    alongside = (along >= t0) & (along <= t1)
    strip = alongside & (dist >= near) & (dist <= far)
    if strip.sum() >= 3:
        return pts[strip]
    return pts[np.argsort(dist)[:20]]


def classify(a: Plane, b: Plane, point, direction, pts_a: np.ndarray, pts_b: np.ndarray) -> str:
    """Ridge/hip: both facets fall away from the edge (convex, water runs off it).
    Valley: both rise away from it (concave, water runs along it).
    Break: anything else, e.g. a steep mansard side meeting a flat top, where only one
    facet slopes. Ridge vs hip: is the line horizontal?"""
    sa = slope_away(a, point, direction, pts_a)
    sb = slope_away(b, point, direction, pts_b)
    if sa < -FLAT_TOL and sb < -FLAT_TOL:
        horizontal = abs(direction[2]) < np.sin(np.radians(HORIZONTAL_TOL_DEG))
        return "ridge" if horizontal else "hip"
    if sa > FLAT_TOL and sb > FLAT_TOL:
        return "valley"
    return "break"


def internal_edges(
    xyz: np.ndarray, seg: Segmentation, radius: float = 0.3, line_band: float = 0.15, snap_dist: float = 1.0
) -> list[Edge]:
    """Every edge where two facets meet.

    1. Adjacent pairs and their contact points (adjacency).
    2. The exact line where the two planes meet (plane_intersection).
    3. The edge's extent: contact points within line_band of that line, projected onto it.
       The band stops points that are near the other facet but not near this edge (for
       example around a hip-end apex) from stretching the edge.
    4. Endpoints: where a third facet touches both, the three planes meet at one point;
       if that corner is within snap_dist of an extent end, use it. It is exact, whereas
       the extent stops roughly half a point spacing short.
    """
    pairs = adjacency(xyz, seg.labels, len(seg.planes), radius)
    edges = []
    for (a, b), contact in pairs.items():
        line = plane_intersection(seg.planes[a], seg.planes[b])
        if line is None:
            continue
        point, d = line
        rel = xyz[contact] - point
        along = rel @ d
        off_line = np.linalg.norm(rel - np.outer(along, d), axis=1)
        on_edge = off_line < line_band
        if on_edge.sum() < 3:
            continue
        t0, t1 = along[on_edge].min(), along[on_edge].max()

        corners = []
        for c in range(len(seg.planes)):
            if c in (a, b) or tuple(sorted((a, c))) not in pairs or tuple(sorted((b, c))) not in pairs:
                continue
            q = three_plane_point(seg.planes[a], seg.planes[b], seg.planes[c])
            if q is not None:
                corners.append((q - point) @ d)
        # Each corner can fix at most one end: when the contact points stop well short, both
        # ends can be within snap_dist of the same corner, which would collapse the edge.
        ends = [t0, t1]
        snapped = [False, False]
        free = list(corners)
        for i in sorted((0, 1), key=lambda i: min((abs(c - ends[i]) for c in free), default=np.inf)):
            if not free:
                break
            best = min(free, key=lambda c: abs(c - ends[i]))
            if abs(best - ends[i]) < snap_dist:
                ends[i], snapped[i] = best, True
                free.remove(best)
        t0, t1 = min(ends), max(ends)
        snapped = tuple(snapped) if ends[0] <= ends[1] else tuple(reversed(snapped))

        pa, pb = (beside_edge(xyz[seg.labels == f], point, d, t0, t1) for f in (a, b))
        kind = classify(seg.planes[a], seg.planes[b], point, d, pa, pb)
        edges.append(Edge(kind, (a, b), point + t0 * d, point + t1 * d, snapped))
    return edges
