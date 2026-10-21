"""Per-plane measurements and a per-house summary of a segmentation."""

from __future__ import annotations

import numpy as np

from roofrecon.planes import azimuth_deg, hull_area, pitch_deg, plan_hull_area
from roofrecon.segment import Segmentation

WALL_PITCH_DEG = 75.0  # steeper than this is a wall, not a roof facet
FLAT_PITCH_DEG = 5.0


def plane_table(xyz: np.ndarray, seg: Segmentation, density: float) -> list[dict]:
    """One row per plane, largest first. density is the house's points per m² of plan area."""
    rows = []
    for k, p in enumerate(seg.planes):
        pts = xyz[seg.labels == k]
        pitch = pitch_deg(p.normal)
        az = azimuth_deg(p.normal)
        plan_area = plan_hull_area(pts)
        rows.append({
            "id": k,
            "points": int(len(pts)),
            "pitch_deg": round(pitch, 2),
            "azimuth_deg": None if az is None else round(az, 1),
            "kind": "wall" if pitch > WALL_PITCH_DEG else ("flat" if pitch < FLAT_PITCH_DEG else "sloped"),
            "rms_m": round(p.rms, 4),
            "true_area_m2": round(hull_area(pts, p), 2),
            "plan_area_m2": round(plan_area, 2),
            # points / (plan area x house density): ~1 for a solid facet. Far below 1 means
            # the points are scattered over a large hull, i.e. not a real facet.
            "fill": round(len(pts) / (plan_area * density), 2) if plan_area > 0 else None,
            "centroid_local": np.round(p.centroid, 3).tolist(),
            "normal": np.round(p.normal, 5).tolist(),
        })
    rows.sort(key=lambda d: -d["points"])
    return rows


def summarize(xyz: np.ndarray, seg: Segmentation, density: float) -> dict:
    planes = plane_table(xyz, seg, density)
    roofs = [d for d in planes if d["kind"] != "wall"]
    return {
        "points": int(len(xyz)),
        "coverage": round(seg.coverage, 4),
        "planes": len(planes),
        "roof_planes": len(roofs),
        "wall_planes": len(planes) - len(roofs),
        "median_rms_m": round(float(np.median([d["rms_m"] for d in planes])), 4) if planes else None,
        "roof_true_area_m2": round(sum(d["true_area_m2"] for d in roofs), 1),
        "per_plane": planes,
    }
