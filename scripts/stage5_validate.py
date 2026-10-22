"""Stage 5: compare each house's report with footprint-based estimates and hand measurements.

    uv run python scripts/stage5_validate.py

Footprint estimates (see roofrecon.validate.footprint_estimates) use the BAG footprint and the
report's pitches, so they check consistency, not correctness. Hand measurements, from
validation/hand_measurements.csv, are the independent check (see validation/PROTOCOL.md).
Writes data/interim/validation.json.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import open3d as o3d
from shapely.geometry import Polygon

from roofrecon.validate import compare_hand_measurements, error_summary, footprint_estimates, overhang, plan_overlap

ROOT = Path("data/interim")
HAND = Path("validation/hand_measurements.csv")
ESTIMATES = [
    "A_footprint_over_cos_mean_pitch",
    "B_footprint_mean_secant",
    "C_widened_footprint_mean_secant",
    "D_footprint_plus_overhang_mean_secant",
]


def main() -> None:
    houses = sorted(p.parent.name for p in ROOT.glob("*/report.json"))
    out = {}
    print(
        f"{'house':10} {'report m²':>9} "
        + " ".join(f"{k.split('_')[0]:>14}" for k in ESTIMATES)
        + "  double-counted m²"
    )
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
        out[h] = {
            "report_area_m2": area,
            "overhang": oh,
            "footprint_estimates": est,
            "plan_overlap": ov,
            "relative_to_report": {k: round(est[k] / area - 1, 4) for k in ESTIMATES},
        }
        print(
            f"{h:10} {area:>9.1f} "
            + " ".join(f"{est[k]:>7.1f} ({100 * (est[k] / area - 1):+5.1f}%)" for k in ESTIMATES)
            + f"  {ov['double_counted_m2']:>5.1f}"
        )
    rows = list(csv.DictReader(HAND.open())) if HAND.exists() else []
    results = []
    for h in sorted({r["house_id"] for r in rows}):
        rep_path = ROOT / h / "report.json"
        if not rep_path.exists():
            print(f"\n{h}: no report.json yet, skipping its hand measurements")
            continue
        results += compare_hand_measurements([r for r in rows if r["house_id"] == h], json.loads(rep_path.read_text()))
    if results:
        print(f"\nHand measurements ({len(results)}):")
        print(f"  {'house':9} {'kind':11} {'feature':34} {'hand':>8} {'report':>8} {'error':>7}  report id/type")
        for r in results:
            rep = "-" if r["report"] is None else f"{r['report']:.2f}"
            err = "-" if r["error"] is None else f"{r['error']:+.2f}"
            print(
                f"  {r['house_id']:9} {r['kind']:11} {r['feature'][:34]:34} {r['hand']:>8.2f} {rep:>8} {err:>7}  "
                f"{r['report_id']} {r.get('report_type') or ''}"
            )
        summary = error_summary(results)
        print("\nError distribution (report minus hand; m for lengths, degrees for slopes):")
        for kind, s_ in summary.items():
            print(
                f"  {kind:11} n={s_['n']:<3} median {s_['median']:+.3f}  "
                f"IQR [{s_['iqr'][0]:+.3f}, {s_['iqr'][1]:+.3f}]  "
                f"max |error| {s_['max_abs']:.3f}"
            )
        repeat = hand_repeatability(results)
        if repeat:
            print(
                "Your own repeatability (|repeat 1 - repeat 2|, same feature): "
                + ", ".join(f"{k} median {v['median']:.3f} max {v['max']:.3f} (n={v['n']})" for k, v in repeat.items())
            )
        out["hand_measurements"] = {"results": results, "summary": summary, "repeatability": repeat}
    (ROOT / "validation.json").write_text(json.dumps(out, indent=2))


def hand_repeatability(results: list[dict]) -> dict:
    """How well hand measurements agree with themselves: same house, kind and feature, repeats 1 and 2."""
    pairs = defaultdict(dict)
    for r in results:
        pairs[(r["house_id"], r["kind"], r["feature"])][str(r.get("repeat"))] = r["hand"]
    diffs = defaultdict(list)
    for (_, kind, _), reps in pairs.items():
        if "1" in reps and "2" in reps:
            diffs[kind].append(abs(reps["1"] - reps["2"]))
    return {
        k: {"n": len(v), "median": round(float(np.median(v)), 3), "max": round(float(max(v)), 3)}
        for k, v in diffs.items()
    }


if __name__ == "__main__":
    main()
