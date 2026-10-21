"""Stage 5: compare each house's report with footprint-based estimates and hand measurements.

    uv run python scripts/stage5_validate.py

Footprint estimates (see roofrecon.validate.footprint_estimates) use the BAG footprint and the
report's pitches, so they check consistency, not correctness. Hand measurements, from
validation/hand_measurements.csv, are the independent check (see validation/PROTOCOL.md).
Writes data/interim/validation.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import open3d as o3d
from shapely.geometry import Polygon

from roofrecon.validate import footprint_estimates, overhang, plan_overlap

ROOT = Path("data/interim")
ESTIMATES = ["A_footprint_over_cos_mean_pitch", "B_footprint_mean_secant", "C_widened_footprint_mean_secant",
             "D_footprint_plus_overhang_mean_secant"]


def main() -> None:
    houses = sorted(p.parent.name for p in ROOT.glob("*/report.json"))
    out = {}
    print(f"{'house':10} {'report m²':>9} " + " ".join(f"{k.split('_')[0]:>14}" for k in ESTIMATES) + "  double-counted m²")
    for h in houses:
        st = json.loads((ROOT / h / "stats.json").read_text())
        rep = json.loads((ROOT / h / "report.json").read_text())
        fp = Polygon(np.array(st["footprint_local_xy"]))
        density = st["density_inside_walls"]["nominal_pts_per_m2"]
        xyz = np.asarray(o3d.io.read_point_cloud(str(ROOT / h / "roof.ply")).points)
        oh = overhang(xyz[:, :2], fp)
        counted = [f for f in rep["facets"] if f["area_m2"] > 0]
        est = footprint_estimates(fp, counted, oh["overhang_m"], oh["points_outside"] / density)
        ov = plan_overlap(counted)
        area = rep["totals"]["roof_area_m2"]
        out[h] = {"report_area_m2": area, "overhang": oh, "footprint_estimates": est, "plan_overlap": ov,
                  "relative_to_report": {k: round(est[k] / area - 1, 4) for k in ESTIMATES}}
        print(f"{h:10} {area:>9.1f} " + " ".join(f"{est[k]:>7.1f} ({100 * (est[k] / area - 1):+5.1f}%)" for k in ESTIMATES)
              + f"  {ov['double_counted_m2']:>5.1f}")
    (ROOT / "validation.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
