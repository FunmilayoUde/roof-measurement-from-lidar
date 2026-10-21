"""Synthetic roofs with known facets and edges, for testing segmentation and topology.

Points are sampled uniformly at random in plan view (like aerial LiDAR: density is per m^2
of ground), then lifted onto the roof surface. Each roof returns per-point facet labels,
the ground-truth plan area of every facet (true area = plan area / cos(pitch)), and the
ground-truth internal edges (where two facets meet), noise-free.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SyntheticRoof:
    xyz: np.ndarray  # (N, 3)
    labels: np.ndarray  # (N,) facet index
    pitch_deg: float
    plan_areas: list[float]  # per facet, m^2 of ground covered
    # Internal edges: {"type": ridge|hip|valley, "facets": (a, b), "p0": xyz, "p1": xyz}
    edges: list[dict] = field(default_factory=list)

    @property
    def true_areas(self) -> list[float]:
        return [a / np.cos(np.radians(self.pitch_deg)) for a in self.plan_areas]


def _sample_plan(rng, length, width, density):
    n = rng.poisson(density * length * width)
    return rng.uniform(0, length, n), rng.uniform(0, width, n)


def gable(length=10.0, width=8.0, pitch_deg=35.0, eave_z=5.0, density=33.0, noise=0.03, seed=0) -> SyntheticRoof:
    """Ridge along x at y = width/2. Facet 0 slopes down to y = 0 (faces south), facet 1 to y = width."""
    rng = np.random.default_rng(seed)
    x, y = _sample_plan(rng, length, width, density)
    labels = (y > width / 2).astype(int)
    dist_to_eave = np.where(labels == 0, y, width - y)
    z = eave_z + np.tan(np.radians(pitch_deg)) * dist_to_eave + rng.normal(0, noise, len(x))
    half = length * width / 2
    ridge_z = eave_z + np.tan(np.radians(pitch_deg)) * width / 2
    edges = [_edge("ridge", (0, 1), (0, width / 2, ridge_z), (length, width / 2, ridge_z))]
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, [half, half], edges)


def hip(length=12.0, width=8.0, pitch_deg=30.0, eave_z=5.0, density=33.0, noise=0.03, seed=0) -> SyntheticRoof:
    """All four sides at the same pitch. With equal pitches, each point belongs to the facet of
    its nearest eave, and the height is eave_z + tan(pitch) * (distance to that eave).
    Facets: 0 south (y = 0), 1 north (y = width), 2 west (x = 0), 3 east (x = length).
    """
    rng = np.random.default_rng(seed)
    x, y = _sample_plan(rng, length, width, density)
    d = np.column_stack([y, width - y, x, length - x])
    labels = d.argmin(axis=1)
    z = eave_z + np.tan(np.radians(pitch_deg)) * d.min(axis=1) + rng.normal(0, noise, len(x))
    ridge = length - width
    trapezoid = (length + ridge) / 2 * (width / 2)
    triangle = width * (width / 2) / 2
    h = width / 2
    zr = eave_z + np.tan(np.radians(pitch_deg)) * h
    a, b = (h, h, zr), (length - h, h, zr)  # ridge ends, where the hips meet
    edges = [
        _edge("ridge", (0, 1), a, b),
        _edge("hip", (0, 2), (0, 0, eave_z), a),
        _edge("hip", (1, 2), (0, width, eave_z), a),
        _edge("hip", (0, 3), (length, 0, eave_z), b),
        _edge("hip", (1, 3), (length, width, eave_z), b),
    ]
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, [trapezoid, trapezoid, triangle, triangle], edges)


def cross_gable(pitch_deg=35.0, eave_z=5.0, density=33.0, noise=0.03, seed=0) -> SyntheticRoof:
    """A T: main gable A over [0,12] x [0,8], ridge along x at y = 4; wing gable B over
    [4,8] x [-6,4], ridge along y at x = 6, same pitch, joining A's south slope.

    Where both roofs cover a spot, the higher one is the visible roof. B's slopes and A's
    south slope meet in two valleys, from the inner corners (4,0) and (8,0) up to (6,2),
    where B's ridge runs into A's slope. A's south facet is notched by B (non-convex).
    Facets: 0 A south, 1 A north, 2 B west, 3 B east.
    """
    rng = np.random.default_rng(seed)
    t = np.tan(np.radians(pitch_deg))
    xa, ya = rng.uniform(0, 12, n := rng.poisson(density * 96)), rng.uniform(0, 8, n)
    xb, yb = rng.uniform(4, 8, n := rng.poisson(density * 24)), rng.uniform(-6, 0, n)
    x, y = np.concatenate([xa, xb]), np.concatenate([ya, yb])
    za = np.where(y <= 4, t * y, t * (8 - y))  # A's height above the eave (-inf outside A)
    za = np.where(y >= 0, za, -np.inf)
    in_b = (x >= 4) & (x <= 8) & (y <= 4)  # B's footprint; without the y bound B would poke through A's north slope
    zb = np.where(in_b, t * np.minimum(x - 4, 8 - x), -np.inf)
    on_b = zb > za
    labels = np.where(on_b, np.where(x < 6, 2, 3), np.where(y <= 4, 0, 1))
    z = eave_z + np.maximum(za, zb) + rng.normal(0, noise, len(x))
    # B covers the part of A where y < min(x-4, 8-x): a triangle of base 4 and height 2.
    b_over_a = 0.5 * 4 * 2
    plan = [48 - b_over_a, 48, (24 + b_over_a) / 2, (24 + b_over_a) / 2]
    e, top = eave_z, eave_z + 2 * t
    edges = [
        _edge("ridge", (0, 1), (0, 4, eave_z + 4 * t), (12, 4, eave_z + 4 * t)),
        _edge("ridge", (2, 3), (6, -6, top), (6, 2, top)),
        _edge("valley", (0, 2), (4, 0, e), (6, 2, top)),
        _edge("valley", (0, 3), (8, 0, e), (6, 2, top)),
    ]
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, plan, edges)


def _edge(kind, facets, p0, p1) -> dict:
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    return {"type": kind, "facets": tuple(facets), "p0": p0, "p1": p1, "length": float(np.linalg.norm(p1 - p0))}


def dormer(pitch_deg=40.0, eave_z=5.0, density=33.0, noise=0.03, seed=0) -> SyntheticRoof:
    """Main gable over [0,12] x [0,9], ridge along x at y = 4.5, with a gable dormer on the
    south slope: 3 m wide (x 4..7), front wall at y = 0.8, eaves at the height the main slope
    reaches at y = 2.0, same pitch, its ridge along y at x = 5.5.

    Seen from above the dormer sits inside the main south slope, so that facet has a hole:
    the dormer's plan area. The dormer's side walls (cheeks) and front wall are vertical and
    get no points. Facets: 0 main south, 1 main north, 2 dormer west, 3 dormer east.
    """
    rng = np.random.default_rng(seed)
    t = np.tan(np.radians(pitch_deg))
    n = rng.poisson(density * 12 * 9)
    x, y = rng.uniform(0, 12, n), rng.uniform(0, 9, n)
    z_main = np.where(y <= 4.5, t * y, t * (9 - y))
    z_d_eave = 2.0 * t
    in_d = (x >= 4) & (x <= 7) & (y >= 0.8) & (y <= 4.5)  # south slope only, or it pokes through the north one
    z_dorm = np.where(in_d, z_d_eave + t * (1.5 - np.abs(x - 5.5)), -np.inf)
    on_d = z_dorm > z_main
    labels = np.where(on_d, np.where(x < 5.5, 2, 3), np.where(y <= 4.5, 0, 1))
    z = eave_z + np.maximum(z_main, z_dorm) + rng.normal(0, noise, n)
    d_facet = 1.2 * 1.5 + 1.5**2 / 2  # x 4..5.5, y from 0.8 up to the valley y = 2 + (x - 4)
    plan = [54.0 - 2 * d_facet, 54.0, d_facet, d_facet]
    e, zd, zdr, zr = eave_z, eave_z + z_d_eave, eave_z + z_d_eave + 1.5 * t, eave_z + 4.5 * t
    edges = [
        _edge("ridge", (0, 1), (0, 4.5, zr), (12, 4.5, zr)),
        _edge("ridge", (2, 3), (5.5, 0.8, zdr), (5.5, 3.5, zdr)),
        _edge("valley", (0, 2), (4, 2.0, zd), (5.5, 3.5, zdr)),
        _edge("valley", (0, 3), (7, 2.0, zd), (5.5, 3.5, zdr)),
    ]
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, plan, edges)
