"""Isolating one building's roof points: class filter, footprint crop, recentering, density."""

from __future__ import annotations

import numpy as np
import shapely
from scipy.spatial import cKDTree
from shapely.geometry import Polygon

# ASPRS LAS 1.4 standard classes (AHN follows these; 26 is AHN's "civil structure").
ASPRS_CLASSES = {
    0: "never classified",
    1: "unclassified",
    2: "ground",
    3: "low vegetation",
    4: "medium vegetation",
    5: "high vegetation",
    6: "building",
    7: "low point (noise)",
    9: "water",
    17: "bridge deck",
    18: "high noise",
    26: "civil structure (AHN)",
}


def classification_histogram(classification: np.ndarray) -> dict[int, int]:
    codes, counts = np.unique(classification, return_counts=True)
    return {int(c): int(n) for c, n in zip(codes, counts)}


def in_polygon(xy: np.ndarray, polygon: Polygon) -> np.ndarray:
    """Boolean mask of 2D points inside polygon (boundary counts as inside)."""
    return shapely.intersects_xy(polygon, xy[:, 0], xy[:, 1])


def crop_mask(xy: np.ndarray, footprint: Polygon, buffer_m: float) -> np.ndarray:
    """Points inside the footprint grown outward by buffer_m.

    Footprints trace the walls; the roof overhangs them by the eave depth. buffer_m has to
    cover that overhang without reaching the neighbour's roof. Mitred joins keep corners
    square instead of rounding them off, which matters because eave corners are square.
    """
    region = footprint.buffer(buffer_m, join_style="mitre") if buffer_m else footprint
    return in_polygon(xy, region)


def buffer_sweep(
    xy: np.ndarray, footprint: Polygon, buffers_m: list[float]
) -> list[dict]:
    """Count building points captured as the crop grows.

    Beyond the eave the count should plateau: no roof, and no building-class points,
    until the next building. A second rise means a neighbour is being swallowed.
    """
    rows = []
    prev = None
    for b in sorted(buffers_m):
        n = int(crop_mask(xy, footprint, b).sum())
        rows.append({"buffer_m": b, "points": n, "added": None if prev is None else n - prev})
        prev = n
    return rows


def recenter(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Shift x, y to a whole-metre local origin; leave z absolute.

    Projected coordinates are large (RD New y ~ 450 000 m, UTM northing ~ 4-5 million m).
    float32, which PLY viewers and GPUs commonly use, has a 24-bit mantissa, so at 450 000 m
    adjacent representable values are ~3 cm apart; at 4.5 million m, ~50 cm. That is
    coarser than the LiDAR itself. Subtracting a whole-metre offset keeps local values
    small (tens of metres) and keeps the offset exactly representable so it can be
    added back losslessly. z is already small (height above datum) and meaningful, so
    it stays put.
    """
    offset = np.array([np.floor(xyz[:, 0].mean()), np.floor(xyz[:, 1].mean()), 0.0])
    return xyz - offset, offset


def density_stats(xy: np.ndarray, area_m2: float, radius_m: float = 1.0) -> dict:
    """Density and spacing, all measured in plan view (per m^2 of ground, the LiDAR-spec
    convention, not per m^2 of sloped roof).

    - nominal: count / area. Depends on what area you divide by, so it inherits any
      error in the footprint or buffer.
    - local: median over points of (neighbours within radius_m) / (pi r^2). Makes no
      assumption about the sampling pattern. Points near the roof edge see half a disc,
      which the median mostly ignores.
    - nearest-neighbour spacing: the resolution you actually have. It does NOT convert
      to density without assuming a pattern (grid: 1/s^2; random: ln2/(pi s^2)), which
      is why it is reported separately.
    """
    n = len(xy)
    out = {"points": n, "area_m2": float(area_m2), "nominal_pts_per_m2": n / area_m2 if area_m2 else None}
    if n >= 2:
        tree = cKDTree(xy)
        d, _ = tree.query(xy, k=2)
        nn = d[:, 1]
        neighbours = tree.query_ball_point(xy, r=radius_m, return_length=True) - 1
        out.update(
            local_pts_per_m2=float(np.median(neighbours) / (np.pi * radius_m**2)),
            local_radius_m=radius_m,
            median_nn_spacing_m=float(np.median(nn)),
            p90_nn_spacing_m=float(np.percentile(nn, 90)),
            duplicate_xy_fraction=float(np.mean(nn == 0)),
        )
    return out
