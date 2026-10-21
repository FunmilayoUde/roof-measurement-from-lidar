"""Stages 2-3 for one house with the chosen settings, so scripts don't each re-wire them."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from roofrecon.polygons import FacetPolygon, roof_polygons
from roofrecon.segment import Segmentation, segment_ransac
from roofrecon.topology import Edge, internal_edges

# Chosen in stage 2 by scripts/stage2_sweep.py (see RESULTS.md).
THRESHOLD_M = 0.12
MIN_AREA_M2 = 1.0


@dataclass
class RoofModel:
    xyz: np.ndarray
    seg: Segmentation
    edges: list[Edge]
    polygons: list[FacetPolygon]
    settings: dict


def build_roof(xyz: np.ndarray, density: float, threshold: float = THRESHOLD_M, min_area: float = MIN_AREA_M2,
               seed: int = 0) -> RoofModel:
    min_points = max(3, int(round(min_area * density)))
    seg = segment_ransac(xyz, threshold=threshold, min_points=min_points, seed=seed)
    edges = internal_edges(xyz, seg)
    polys = roof_polygons(xyz, seg, edges)
    settings = {"method": "ransac", "threshold_m": threshold, "min_area_m2": min_area, "min_points": min_points,
                "density_pts_per_m2": round(density, 2), "seed": seed}
    return RoofModel(xyz, seg, edges, polys, settings)
