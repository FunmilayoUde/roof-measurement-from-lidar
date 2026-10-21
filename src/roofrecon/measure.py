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
                edges.append({"type": "unclassified", "facets": [P.facet], "p0": s.p0, "p1": s.p1, "disagreement_m": None})
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


def measurement_report(model: RoofModel, meta: dict) -> dict:
    facets = []
    for P in model.polygons:
        pitch = pitch_deg(P.plane.normal)
        az = azimuth_deg(P.plane.normal)
        facets.append({
            "id": P.facet,
            "points": int((model.seg.labels == P.facet).sum()),
            "pitch_deg": round(pitch, 1),
            "pitch_rise_per_12": round(rise_per_12(pitch), 1),
            "azimuth_deg": None if az is None else round(az, 1),
            "facing": compass(az),
            "area_m2": round(P.area, 2),
            "area_ft2": round(P.area * M2_TO_FT2, 1),
            "plan_area_m2": round(P.plan_area, 2),
            "outline": "traced (fallback)" if P.fallback else "plane and line intersections",
            "fallback_reason": P.fallback_reason,
            "vertices_local_m": np.round(P.vertices, 3).tolist(),
        })

    edges = []
    for i, e in enumerate(sorted(reconcile_edges(model), key=lambda e: (EDGE_TYPES.index(e["type"]), e["facets"]))):
        length = float(np.linalg.norm(e["p1"] - e["p0"]))
        edges.append({
            "id": i, "type": e["type"], "facets": [int(f) for f in e["facets"]],
            "length_m": round(length, 2), "length_ft": round(length * M_TO_FT, 1),
            "endpoint_disagreement_m": None if e["disagreement_m"] is None else round(float(e["disagreement_m"]), 2),
            "p0_local_m": np.round(e["p0"], 3).tolist(), "p1_local_m": np.round(e["p1"], 3).tolist(),
        })

    area = sum(f["area_m2"] for f in facets)
    by_type = defaultdict(float)
    for e in edges:
        by_type[e["type"]] += e["length_m"]
    fallback = [f["id"] for f in facets if f["outline"] != "plane and line intersections"]
    notes = []
    if fallback:
        notes.append(f"{len(fallback)} facet(s) use a traced outline, not line intersections: their area is "
                     "approximate and their sides are unclassified.")
    disagreements = [e["endpoint_disagreement_m"] for e in edges if e["endpoint_disagreement_m"] is not None]
    if disagreements and max(disagreements) > 0.25:
        notes.append(f"Shared edges were reconciled from two facets that disagreed by up to {max(disagreements):.2f} m "
                     "at an endpoint (above 0.5 m, the overlap of the two versions is kept).")
    return {
        "house_id": meta.get("house_id"),
        "source": meta,
        "settings": model.settings,
        "units": {"length": "m (and ft)", "area": "m² (and ft²)", "square": "100 ft²"},
        "totals": {
            "roof_area_m2": round(area, 1),
            "roof_area_ft2": round(area * M2_TO_FT2, 0),
            "squares": round(area * M2_TO_FT2 / FT2_PER_SQUARE, 2),
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
