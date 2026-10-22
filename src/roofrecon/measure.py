"""Stage 4: the measurement report for one roof.

Facets come from their polygons. Edges are listed once each: a shared internal edge appears
as a side of both facet polygons, usually with slightly different endpoints (each facet ends
it where its own outer lines cross it), so the two versions are reconciled by averaging
matching endpoints.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from roofrecon.pipeline import RoofModel
from roofrecon.planes import azimuth_deg, pitch_deg

M_TO_FT = 3.280839895
M2_TO_FT2 = M_TO_FT**2
FT2_PER_SQUARE = 100.0
COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
EDGE_TYPES = ["ridge", "hip", "valley", "break", "eave", "rake", "top", "perimeter", "unclassified"]

# Facet area error by size, measured on the synthetic roofs (gable, hip, cross-gable, dormer;
# 20 seeds each; RESULTS.md, stage 5): 90th percentile of |area error|. Those roofs have no
# facets between 3.8 and 17.1 m², so the 10 m² boundary could sit anywhere in that range.
# Synthetic roofs are ideal (exact planes, no trees), so real errors will be at least this.
SMALL_FACET_M2 = 10.0
AREA_ERR_P90 = [(SMALL_FACET_M2, 0.315), (25.0, 0.012), (float("inf"), 0.007)]


def area_uncertainty(area_m2: float, fallback: bool) -> tuple[str, float]:
    """(confidence, 90th-percentile relative area error) for a facet of this size."""
    p90 = next(e for limit, e in AREA_ERR_P90 if area_m2 < limit)
    low = fallback or area_m2 < SMALL_FACET_M2
    return ("low" if low else "normal"), (max(p90, AREA_ERR_P90[0][1]) if fallback else p90)


def compass(az: float | None) -> str | None:
    return None if az is None else COMPASS[int(((az + 11.25) % 360) // 22.5)]


def rise_per_12(pitch: float) -> float:
    """US roofing pitch: inches of rise per 12 inches of run."""
    return 12.0 * np.tan(np.radians(pitch))


def _endpoint_pairs(a0, a1, b0, b1, overlap_above: float = 0.5):
    """Combine two facets' versions of one shared edge (same line, possibly different extents).

    Small disagreement is noise in where each facet's outer lines cross the edge: average the
    matched endpoints. Large disagreement means the versions cover different stretches, e.g.
    one facet's ridge also runs along a third facet (delft_05: 3.5 m); then keep the overlap,
    the stretch where both facets actually border each other. Returns (p0, p1, disagreement).
    """
    same = np.linalg.norm(a0 - b0) + np.linalg.norm(a1 - b1)
    swap = np.linalg.norm(a0 - b1) + np.linalg.norm(a1 - b0)
    if swap < same:
        b0, b1 = b1, b0
    gap = max(np.linalg.norm(a0 - b0), np.linalg.norm(a1 - b1))
    if gap <= overlap_above:
        return (a0 + b0) / 2, (a1 + b1) / 2, gap
    d = (a1 - a0) / np.linalg.norm(a1 - a0)
    ta, tb = (0.0, (a1 - a0) @ d), sorted(((b0 - a0) @ d, (b1 - a0) @ d))
    lo, hi = max(ta[0], tb[0]), min(ta[1], tb[1])
    if hi <= lo:  # no overlap at all: fall back to the average
        return (a0 + b0) / 2, (a1 + b1) / 2, gap
    return a0 + lo * d, a0 + hi * d, gap


def reconcile_edges(model: RoofModel) -> list[dict]:
    """One record per physical edge.

    Internal sides are grouped by facet pair and type; when both facets have the side, the two
    versions are combined (_endpoint_pairs) and the disagreement reported. Outer sides belong
    to one facet. Sides of fallback facets (traced outline) and outer sides that are neither
    level nor up-slope are 'unclassified'.
    """
    internal = defaultdict(list)
    edges = []
    for P in model.polygons:
        for s in P.sides:
            if P.fallback:
                edges.append(
                    {"type": "unclassified", "facets": [P.facet], "p0": s.p0, "p1": s.p1, "disagreement_m": None}
                )
            elif s.neighbour is None:
                kind = "unclassified" if s.kind == "other" else s.kind  # outer side neither level nor up-slope
                edges.append({"type": kind, "facets": [P.facet], "p0": s.p0, "p1": s.p1, "disagreement_m": None})
            else:
                internal[(min(P.facet, s.neighbour), max(P.facet, s.neighbour), s.kind)].append(s)
    for (a, b, kind), sides in internal.items():
        if len(sides) >= 2:
            p0, p1, gap = _endpoint_pairs(sides[0].p0, sides[0].p1, sides[1].p0, sides[1].p1)
        else:  # the other facet fell back, or didn't produce this side
            p0, p1, gap = sides[0].p0, sides[0].p1, None
        edges.append({"type": kind, "facets": [a, b], "p0": p0, "p1": p1, "disagreement_m": gap})
    return edges


def trim_fallback_overlaps(model: RoofModel, max_dz: float = 0.5) -> dict[int, float]:
    """Area of each fallback facet that should not be counted again.

    A fallback facet (traced outline) can lie inside a facet closed by line intersections at
    the same height: on delft_01 two traced pieces of the curved top (24 m²) sat inside the flat
    top's polygon and were counted twice. Overlap with a facet whose surface is more than
    max_dz higher or lower there is kept: that's a higher roof overhanging a lower one, and
    the lower roof really does continue underneath.
    Returns {facet: true (sloped) area to subtract}.
    """
    trims = {}
    good = [P for P in model.polygons if not P.fallback]
    for F in (P for P in model.polygons if P.fallback):
        f_plan = F.plan_shape()
        overlap = 0.0
        for G in good:
            inter = f_plan.intersection(G.plan_shape())
            if inter.area < 1e-6:
                continue
            cx, cy = inter.centroid.x, inter.centroid.y
            z_f, z_g = (_height_at(P.plane, cx, cy) for P in (F, G))
            if abs(z_f - z_g) <= max_dz:
                overlap += inter.area
        if overlap > 0:
            frac = min(1.0, overlap / max(f_plan.area, 1e-9))
            trims[F.facet] = frac * F.area
    return trims


def _height_at(plane, x: float, y: float) -> float:
    n, c = plane.normal, plane.centroid
    return float(c[2] - (n[0] * (x - c[0]) + n[1] * (y - c[1])) / n[2]) if abs(n[2]) > 1e-6 else float(c[2])


def measurement_report(model: RoofModel, meta: dict) -> dict:
    trims = trim_fallback_overlaps(model)
    facets = []
    for P in model.polygons:
        pitch = pitch_deg(P.plane.normal)
        az = azimuth_deg(P.plane.normal)
        area = P.area - trims.get(P.facet, 0.0)
        confidence, err_p90 = area_uncertainty(area, P.fallback)
        facets.append(
            {
                "id": P.facet,
                "points": int((model.seg.labels == P.facet).sum()),
                "pitch_deg": round(pitch, 1),
                "pitch_rise_per_12": round(rise_per_12(pitch), 1),
                "azimuth_deg": None if az is None else round(az, 1),
                "facing": compass(az),
                "area_m2": round(area, 2),
                "area_ft2": round(area * M2_TO_FT2, 1),
                "area_not_counted_m2": round(trims[P.facet], 2) if P.facet in trims else None,
                "confidence": confidence,
                "area_error_p90_pct": round(100 * err_p90, 1),
                "plan_area_m2": round(P.plan_area, 2),
                "outline": "traced (fallback)" if P.fallback else "plane and line intersections",
                "fallback_reason": P.fallback_reason,
                "vertices_local_m": np.round(P.vertices, 3).tolist(),
                "holes_local_m": [np.round(h, 3).tolist() for h in P.holes],
            }
        )

    edges = []
    for i, e in enumerate(sorted(reconcile_edges(model), key=lambda e: (EDGE_TYPES.index(e["type"]), e["facets"]))):
        length = float(np.linalg.norm(e["p1"] - e["p0"]))
        edges.append(
            {
                "id": i,
                "type": e["type"],
                "facets": [int(f) for f in e["facets"]],
                "length_m": round(length, 2),
                "length_ft": round(length * M_TO_FT, 1),
                "endpoint_disagreement_m": None
                if e["disagreement_m"] is None
                else round(float(e["disagreement_m"]), 2),
                "p0_local_m": np.round(e["p0"], 3).tolist(),
                "p1_local_m": np.round(e["p1"], 3).tolist(),
            }
        )

    area = sum(f["area_m2"] for f in facets)
    low_area = sum(f["area_m2"] for f in facets if f["confidence"] == "low" and f["area_m2"] > 0)
    # Conservative: facet errors added as if they all went the same way.
    area_bound = sum(f["area_m2"] * f["area_error_p90_pct"] / 100 for f in facets)
    by_type = defaultdict(float)
    for e in edges:
        by_type[e["type"]] += e["length_m"]
    fallback = [f["id"] for f in facets if f["outline"] != "plane and line intersections"]
    notes = []
    if fallback:
        notes.append(
            f"{len(fallback)} facet(s) use a traced outline, not line intersections: their area is "
            "approximate and their sides are unclassified."
        )
    if low_area > 0:
        notes.append(
            f"{low_area:.1f} m² ({100 * low_area / area:.0f}% of the roof) is in small or traced facets, whose "
            f"area error was up to {100 * AREA_ERR_P90[0][1]:.0f}% (90th percentile) on synthetic roofs; "
            "real errors will be at least that."
        )
    if trims:
        notes.append(
            f"{sum(trims.values()):.1f} m² of traced facet(s) {sorted(trims)} lies inside other facets at the "
            "same height and is not counted twice."
        )
    disagreements = [e["endpoint_disagreement_m"] for e in edges if e["endpoint_disagreement_m"] is not None]
    if disagreements and max(disagreements) > 0.25:
        notes.append(
            f"Shared edges were reconciled from two facets that disagreed by up to {max(disagreements):.2f} m "
            "at an endpoint (above 0.5 m, the overlap of the two versions is kept)."
        )
    return {
        "house_id": meta.get("house_id"),
        "source": meta,
        "settings": model.settings,
        "units": {"length": "m (and ft)", "area": "m² (and ft²)", "square": "100 ft²"},
        "totals": {
            "roof_area_m2": round(area, 1),
            "roof_area_ft2": round(area * M2_TO_FT2, 0),
            "squares": round(area * M2_TO_FT2 / FT2_PER_SQUARE, 2),
            "roof_area_p90_error_m2": round(area_bound, 1),
            "area_in_low_confidence_facets_m2": round(low_area, 1),
            "facets": len(facets),
            "length_m_by_type": {t: round(by_type[t], 2) for t in EDGE_TYPES if t in by_type},
            "length_ft_by_type": {t: round(by_type[t] * M_TO_FT, 1) for t in EDGE_TYPES if t in by_type},
        },
        "quality": {
            "coverage": round(model.seg.coverage, 4),
            "unassigned_points": int((model.seg.labels < 0).sum()),
            "fallback_facets": fallback,
            "notes": notes,
        },
        "facets": facets,
        "edges": edges,
    }
