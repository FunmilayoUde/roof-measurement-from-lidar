"""Per-point normals from local neighbourhoods."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def estimate_normals(xyz: np.ndarray, radius: float, max_nn: int = 64) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Normal, surface variation and neighbour count for every point, by local PCA.

    For each point, take its neighbours within `radius` and fit a plane to them (the same
    SVD idea as fit_plane, done here via the 3x3 covariance's eigenvectors). The direction
    of least spread is the normal.

    Surface variation = l0 / (l0 + l1 + l2), with l0 the smallest eigenvalue: 0 on a perfect
    plane, up to 1/3 for points scattered equally in all directions. It is high where the
    neighbourhood straddles a ridge or edge, which is exactly where the normal is
    unreliable.

    The radius trades noise against blur: with ~0.09 m spacing, r = 0.15 m gives only 2-3
    neighbours (normals dominated by noise), while large radii smear normals across ridges
    over a band about r wide. Points with fewer than 3 neighbours get a NaN normal.
    """
    n = len(xyz)
    k = min(max_nn, n)
    dist, idx = cKDTree(xyz).query(xyz, k=k)
    within = dist <= radius  # (n, k): neighbours that are inside the radius
    count = within.sum(axis=1)

    nbrs = xyz[idx]  # (n, k, 3)
    w = within[..., None]
    mean = (nbrs * w).sum(axis=1) / np.maximum(count, 1)[:, None]
    centred = (nbrs - mean[:, None, :]) * w
    cov = np.einsum("nki,nkj->nij", centred, centred) / np.maximum(count, 1)[:, None, None]
    evals, evecs = np.linalg.eigh(cov)  # ascending eigenvalues
    normals = evecs[:, :, 0]
    normals[normals[:, 2] < 0] *= -1
    variation = evals[:, 0] / np.maximum(evals.sum(axis=1), 1e-12)

    bad = count < 3
    normals[bad] = np.nan
    variation[bad] = np.nan
    return normals, variation, count
