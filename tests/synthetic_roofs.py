"""Synthetic roofs with known facets, for testing plane fitting and segmentation.

Points are sampled uniformly at random in plan view (like aerial LiDAR: density is per m^2
of ground), then lifted onto the roof surface. Each roof returns per-point facet labels and
the ground-truth plan area of every facet, so true area = plan area / cos(pitch).
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class SyntheticRoof:
    xyz: np.ndarray  # (N, 3)
    labels: np.ndarray  # (N,) facet index
    pitch_deg: float
    plan_areas: list[float]  # per facet, m^2 of ground covered

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
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, [half, half])


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
    return SyntheticRoof(np.column_stack([x, y, z]), labels, pitch_deg, [trapezoid, trapezoid, triangle, triangle])
