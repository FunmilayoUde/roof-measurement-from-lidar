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


def facet_plan_shape(f: dict):
    """A report facet seen from above, holes removed."""
    outer = shapely.make_valid(Polygon(np.asarray(f["vertices_local_m"])[:, :2]))
    holes = [shapely.make_valid(Polygon(np.asarray(h)[:, :2])) for h in f.get("holes_local_m", [])]
    return outer.difference(shapely.union_all(holes)) if holes else outer


def plan_overlap(facets: list[dict]) -> dict:
    """Double counting: the sum of facet plan areas minus the area of their union. Facets of
    one roof can't legitimately overlap in plan, except where a higher roof overhangs a lower
    one, so anything above that is area counted twice."""
    polys = [facet_plan_shape(f) for f in facets]
    total = sum(p.area for p in polys)
    union = shapely.union_all(polys).area
    pairs = []
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            a = polys[i].intersection(polys[j]).area
            if a > 0.5:
                pairs.append({"facets": [facets[i]["id"], facets[j]["id"]], "overlap_m2": round(a, 1)})
    # Area the report already declined to count twice (trimmed fallback overlaps), in plan.
    trimmed = sum((f.get("area_not_counted_m2") or 0.0) * np.cos(np.radians(f["pitch_deg"])) for f in facets)
    return {"plan_sum_m2": round(total, 1), "plan_union_m2": round(union, 1),
            "already_trimmed_m2": round(trimmed, 1), "double_counted_m2": round(total - union - trimmed, 1),
            "overlapping_pairs": sorted(pairs, key=lambda p: -p["overlap_m2"])}


def _match_edge(p0, p1, edges: list[dict], max_angle_deg: float = 20.0):
    """The report edge nearest a picked segment: smallest midpoint distance among edges within
    max_angle_deg of its direction."""
    mid, d = (p0 + p1) / 2, (p1 - p0) / np.linalg.norm(p1 - p0)
    best, best_dist = None, np.inf
    for e in edges:
        a, b = np.asarray(e["p0_local_m"]), np.asarray(e["p1_local_m"])
        ed = (b - a) / max(np.linalg.norm(b - a), 1e-9)
        if abs(d @ ed) < np.cos(np.radians(max_angle_deg)):
            continue
        dist = np.linalg.norm(mid - (a + b) / 2)
        if dist < best_dist:
            best, best_dist = e, dist
    return best, best_dist


def _match_facet(p0, p1, facets: list[dict]):
    """The facet whose plan polygon contains the picked pair's midpoint (nearest if none does)."""
    mid = (p0 + p1) / 2
    pt = shapely.points(mid[:2])
    scored = []
    for f in facets:
        if f["area_m2"] <= 0:
            continue
        poly = facet_plan_shape(f)
        scored.append((0.0 if poly.contains(pt) else poly.distance(pt), f))
    return min(scored, key=lambda s: s[0])[1] if scored else None


def compare_hand_measurements(rows: list[dict], report: dict) -> list[dict]:
    """One result per measurement row: the hand value, the report's value, and the error
    (report minus hand). Lengths in m, pitch in degrees."""
    edges = {e["id"]: e for e in report["edges"]}
    facets = {f["id"]: f for f in report["facets"]}
    out = []
    for r in rows:
        kind = r["kind"]
        rid = int(r["report_id"]) if r.get("report_id") not in (None, "") else None
        res = {"house_id": r["house_id"], "kind": kind, "feature": r["feature"], "method": r.get("method"),
               "repeat": r.get("repeat"), "report_id": rid}
        if kind in ("edge", "slope"):
            p0 = np.array([float(r[k]) for k in ("x0", "y0", "z0")])
            p1 = np.array([float(r[k]) for k in ("x1", "y1", "z1")])
        if kind == "edge":
            hand = float(np.linalg.norm(p1 - p0))
            e = edges.get(rid) if rid is not None else _match_edge(p0, p1, report["edges"])[0]
            res.update(hand=round(hand, 3), report=None if e is None else e["length_m"],
                       report_id=None if e is None else e["id"], report_type=None if e is None else e["type"])
        elif kind == "slope":
            dz = abs(p1[2] - p0[2])
            run = float(np.linalg.norm((p1 - p0)[:2]))
            hand = float(np.degrees(np.arctan2(dz, run)))
            f = facets.get(rid) if rid is not None else _match_facet(p0, p1, report["facets"])
            res.update(hand=round(hand, 2), report=None if f is None else f["pitch_deg"],
                       report_id=None if f is None else f["id"], report_type="pitch")
        elif kind == "plan_length":
            hand = float(r["value"])
            e = edges.get(rid)
            plan = None if e is None else float(np.linalg.norm((np.asarray(e["p1_local_m"]) - np.asarray(e["p0_local_m"]))[:2]))
            res.update(hand=round(hand, 3), report=None if plan is None else round(plan, 3),
                       report_type=None if e is None else e["type"] + " (plan)")
        else:
            raise ValueError(f"unknown kind {kind!r}")
        res["error"] = None if res["report"] is None else round(res["report"] - res["hand"], 3)
        out.append(res)
    return out


def error_summary(results: list[dict]) -> dict:
    """Distribution of errors per kind, not just a mean: n, median, IQR, max |error|."""
    summary = {}
    for kind in sorted({r["kind"] for r in results}):
        e = np.array([r["error"] for r in results if r["kind"] == kind and r["error"] is not None])
        if len(e) == 0:
            continue
        q1, med, q3 = np.percentile(e, [25, 50, 75])
        summary[kind] = {"n": int(len(e)), "median": round(float(med), 3), "iqr": [round(float(q1), 3), round(float(q3), 3)],
                         "max_abs": round(float(np.abs(e).max()), 3), "mean_abs": round(float(np.abs(e).mean()), 3)}
    return summary
