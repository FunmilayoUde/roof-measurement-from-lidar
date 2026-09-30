"""A synthetic scene with known answers, placed at real RD New coordinates so precision bugs show up."""

import numpy as np
import pytest
from shapely.geometry import box

ORIGIN = np.array([85_000.0, 446_000.0])  # roughly Delft in EPSG:28992
WALLS = (10.0, 8.0)  # footprint x, y extent in metres
OVERHANG = 0.4
EAVE_Z = 5.0
PITCH_DEG = 35.0
NEIGHBOUR_GAP = 3.0  # wall-to-wall gap to the next house (its eaves overhang too)
DENSITY = 20.0  # points per m^2 of plan area


def gable(rng, x0, y0, lx, ly, density):
    """Gable roof with ridge along x, sampled uniformly at random in plan (Poisson)."""
    xmin, ymin = x0 - OVERHANG, y0 - OVERHANG
    xmax, ymax = x0 + lx + OVERHANG, y0 + ly + OVERHANG
    n = rng.poisson(density * (xmax - xmin) * (ymax - ymin))
    x = rng.uniform(xmin, xmax, n)
    y = rng.uniform(ymin, ymax, n)
    ridge_y = y0 + ly / 2
    # Height rises from the wall line at EAVE_Z; the overhang continues the same plane downward.
    z = EAVE_Z + np.tan(np.radians(PITCH_DEG)) * (ly / 2 - np.abs(y - ridge_y))
    return np.column_stack([x, y, z])


@pytest.fixture
def scene():
    rng = np.random.default_rng(0)
    x0, y0 = ORIGIN
    house = gable(rng, x0, y0, *WALLS, DENSITY)
    neighbour = gable(rng, x0 + WALLS[0] + NEIGHBOUR_GAP, y0, *WALLS, DENSITY)

    n_ground = 4000
    ground = np.column_stack([
        rng.uniform(x0 - 5, x0 + 30, n_ground),
        rng.uniform(y0 - 5, y0 + 15, n_ground),
        np.zeros(n_ground),
    ])
    tree = np.column_stack([rng.normal(x0 + 1, 0.5, 300), rng.normal(y0 - 2, 0.5, 300), rng.uniform(2, 8, 300)])

    xyz = np.vstack([house, neighbour, ground, tree])
    cls = np.concatenate([
        np.full(len(house), 6), np.full(len(neighbour), 6), np.full(n_ground, 2), np.full(len(tree), 5)
    ]).astype(np.uint8)
    return {
        "xyz": xyz,
        "classification": cls,
        "n_house": len(house),
        "footprint": box(x0, y0, x0 + WALLS[0], y0 + WALLS[1]),
    }
