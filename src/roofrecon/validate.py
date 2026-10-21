"""Stage 5: checks of the measurement report against independent estimates."""

from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import Polygon


def overhang(xy: np.ndarray, footprint: Polygon, q: float = 98.0) -> dict:
    """How far the roof extends past the wall outline: a high percentile of the distance from
    the outline of roof points that lie outside it. A percentile, not the max, so a few stray
    points don't set it."""
    outside = ~shapely.contains_xy(footprint, xy[:, 0], xy[:, 1])
    d = shapely.distance(footprint, shapely.points(xy[outside]))
    return {"overhang_m": float(np.percentile(d, q)) if len(d) else 0.0, "points_outside": int(outside.sum()),
            "median_outside_m": float(np.median(d)) if len(d) else 0.0}


def footprint_estimates(footprint: Polygon, facets: list[dict], overhang_m: float, overhang_area_m2: float) -> dict:
    """Roof area from the building footprint and the facets' pitches, three ways.

    A: footprint / cos(mean pitch). Mean pitch weighted by plan area.
    B: footprint x mean of 1/cos(pitch), same weights. Correct for a roof of mixed pitches,
       since 1/cos is convex: 1/cos(mean pitch) < mean of 1/cos(pitch) (Jensen).
    C: B with the footprint widened by the overhang (mitred), because the footprint traces
       the walls and the roof extends past them. Assumes the same overhang all round, which
       roofs rarely have (eaves overhang, gable ends and flat roofs barely), so it overshoots.
    D: B with the footprint plus the overhang's plan area, measured as (roof points outside the
       footprint) / density.
    All use the facets' pitches and plan-area shares from the report, but not its areas.
    """
    plan = np.array([f["plan_area_m2"] for f in facets])
    pitch = np.radians([f["pitch_deg"] for f in facets])
    w = plan / plan.sum()
    mean_pitch = float(w @ pitch)
    mean_secant = float(w @ (1 / np.cos(pitch)))
    widened = footprint.buffer(overhang_m, join_style="mitre").area
    return {
        "mean_pitch_deg": round(np.degrees(mean_pitch), 1),
        "secant_of_mean_pitch": round(1 / np.cos(mean_pitch), 4),
        "mean_secant": round(mean_secant, 4),
        "footprint_m2": round(footprint.area, 1),
        "footprint_widened_m2": round(widened, 1),
        "A_footprint_over_cos_mean_pitch": round(footprint.area / np.cos(mean_pitch), 1),
        "B_footprint_mean_secant": round(footprint.area * mean_secant, 1),
        "C_widened_footprint_mean_secant": round(widened * mean_secant, 1),
        "overhang_area_m2": round(overhang_area_m2, 1),
        "D_footprint_plus_overhang_mean_secant": round((footprint.area + overhang_area_m2) * mean_secant, 1),
        "facets_plan_area_sum_m2": round(float(plan.sum()), 1),
    }


def plan_overlap(facets: list[dict]) -> dict:
    """Double counting: the sum of facet plan areas minus the area of their union. Facets of
    one roof can't legitimately overlap in plan, except where a higher roof overhangs a lower
    one, so anything above that is area counted twice."""
    polys = [shapely.make_valid(Polygon(np.asarray(f["vertices_local_m"])[:, :2])) for f in facets]
    total = sum(p.area for p in polys)
    union = shapely.union_all(polys).area
    pairs = []
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            a = polys[i].intersection(polys[j]).area
            if a > 0.5:
                pairs.append({"facets": [facets[i]["id"], facets[j]["id"]], "overlap_m2": round(a, 1)})
    return {"plan_sum_m2": round(total, 1), "plan_union_m2": round(union, 1), "double_counted_m2": round(total - union, 1),
            "overlapping_pairs": sorted(pairs, key=lambda p: -p["overlap_m2"])}
