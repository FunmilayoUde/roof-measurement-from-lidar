"""Plane geometry: fitting, pitch, facing direction, and true surface area."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import shapely


@dataclass
class Plane:
    normal: np.ndarray  # unit vector, pointing up (nz >= 0)
    centroid: np.ndarray  # a point on the plane
    rms: float  # root-mean-square point-to-plane distance of the fitted points, metres

    @property
    def d(self) -> float:
        """Offset in the plane equation n . p + d = 0."""
        return -float(self.normal @ self.centroid)

    def distance(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance of points from the plane (positive = above, since the normal points up)."""
        return (xyz - self.centroid) @ self.normal


def fit_plane(xyz: np.ndarray) -> Plane:
    """Least-squares plane through points, by SVD.

    Centre the points; the plane's normal is the direction in which they spread least,
    which is the singular vector with the smallest singular value. This minimises the
    perpendicular distances to the plane, unlike fitting z = ax + by + c, which minimises
    vertical distances and is biased on steep facets.
    """
    if len(xyz) < 3:
        raise ValueError("Need at least 3 points to fit a plane")
    centroid = xyz.mean(axis=0)
    _, s, vt = np.linalg.svd(xyz - centroid, full_matrices=False)
    normal = vt[-1]
    if normal[2] < 0:  # a normal and its negative describe the same plane; pick the upward one
        normal = -normal
    rms = float(s[-1] / np.sqrt(len(xyz)))  # smallest singular value = residual spread along the normal
    return Plane(normal, centroid, rms)


def pitch_deg(normal: np.ndarray) -> float:
    """Slope of the plane from horizontal: the angle between its normal and vertical."""
    return float(np.degrees(np.arccos(np.clip(abs(normal[2]), 0.0, 1.0))))


def azimuth_deg(normal: np.ndarray) -> float | None:
    """Compass direction the facet faces (downhill), 0 = north (+y), 90 = east (+x).

    The horizontal part of an upward normal points downhill. Undefined for flat planes.
    """
    nx, ny = normal[0], normal[1]
    if np.hypot(nx, ny) < 1e-6:
        return None
    return float(np.degrees(np.arctan2(nx, ny)) % 360)


def plane_axes(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two unit vectors spanning the plane: u is horizontal (along the eave/ridge direction),
    v points up the slope. For a flat plane, u = x and v = y.
    """
    horizontal = np.array([normal[0], normal[1], 0.0])
    if np.linalg.norm(horizontal) < 1e-6:
        return np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    u = np.cross([0.0, 0.0, 1.0], normal)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    return u, v


def to_plane_coords(xyz: np.ndarray, plane: Plane) -> np.ndarray:
    """2D coordinates of points within the plane (after dropping their distance from it)."""
    u, v = plane_axes(plane.normal)
    rel = xyz - plane.centroid
    return np.column_stack([rel @ u, rel @ v])


def hull_area(xyz: np.ndarray, plane: Plane) -> float:
    """True (sloped) surface area: area of the convex hull of the points, measured in the plane.

    Two known biases:
    - Too small: the outermost points sit a little inside the real edge (on average about
      half a point spacing), so the hull is slightly smaller than the facet.
    - Too big: for non-convex facets (an L shape, a notch around a dormer) the hull bridges
      across the gap. Stage 3 replaces this with the facet polygon from plane intersections.
    """
    uv = to_plane_coords(xyz, plane)
    return float(shapely.convex_hull(shapely.multipoints(uv)).area)


def plan_hull_area(xyz: np.ndarray) -> float:
    """Area seen from above (convex hull of x, y). For a planar facet: hull_area * cos(pitch)."""
    return float(shapely.convex_hull(shapely.multipoints(xyz[:, :2])).area)


def alpha_shape(uv: np.ndarray, radius: float):
    """Outline of 2D points: the union of Delaunay triangles whose circumradius is <= radius.

    Unlike a convex hull it follows notches; unlike shapely's concave_hull, whose setting
    is relative to the longest edges, `radius` is an absolute length, so it means the same
    thing on every facet. Too small (near the point spacing) and it bites holes into the
    facet; too large and it bridges real notches. Returns a shapely (Multi)Polygon.
    """
    from scipy.spatial import Delaunay

    tri = Delaunay(uv)
    a, b, c = (uv[tri.simplices[:, i]] for i in range(3))
    la, lb, lc = (np.linalg.norm(p - q, axis=1) for p, q in ((b, c), (a, c), (a, b)))
    area2 = np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    with np.errstate(divide="ignore", invalid="ignore"):
        circumradius = la * lb * lc / (2 * area2)
    keep = tri.simplices[circumradius <= radius]
    return shapely.union_all(shapely.polygons(uv[keep]))
