"""Collect every number and figure input of the project write-up into one JSON file.

    uv run python scripts/make_report_data.py

Reads the stage outputs under data/interim (points, report.json, validation.json, sweeps) and
re-runs the synthetic area check, then writes docs/report/report_data.json. The page
docs/report/index.html draws all its figures and its 3D view from that file, so after new
houses, new settings or hand measurements, re-running this script updates the write-up.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import open3d as o3d

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from roofrecon.export import roof_mesh  # noqa: E402
from roofrecon.measure import measurement_report  # noqa: E402
from roofrecon.pipeline import build_roof  # noqa: E402
from synthetic_roofs import cross_gable, dormer, gable, hip  # noqa: E402

ROOT = Path("data/interim")
OUT = Path("docs/report/report_data.json")
HOUSES = ["delft_03", "delft_02", "delft_05", "delft_01"]
# From the BAG record and the stage 1 preview (RESULTS.md, stage 1).
ABOUT = {
    "delft_01": {"address": "(address withheld)", "built": "1910s", "roof": "mansard with a curved top, flat annex"},
    "delft_02": {"address": "(address withheld)", "built": "1970s", "roof": "two flat levels"},
    "delft_03": {"address": "(address withheld)", "built": "1920s", "roof": "gable, flat extension, small lean-to"},
    "delft_05": {"address": "(address withheld)", "built": "2010s", "roof": "hipped, with notched corners"},
}


def r(a, nd=3):
    return np.round(np.asarray(a, float), nd).tolist()


def house(h: str) -> dict:
    d = ROOT / h
    st = json.loads((d / "stats.json").read_text())
    rep = json.loads((d / "report.json").read_text())
    xyz = np.asarray(o3d.io.read_point_cloud(str(d / "roof.ply")).points)
    density = st["density_inside_walls"]["nominal_pts_per_m2"]
    model = build_roof(xyz, density, seed=rep["settings"]["seed"])
    mesh = roof_mesh(rep)
    return {
        "id": h, **ABOUT[h], "bag_id": rep["source"]["bag_id"],
        "points": len(xyz), "density_pts_per_m2": round(density, 1),
        "nn_spacing_m": round(st["density_inside_walls"]["median_nn_spacing_m"], 3),
        "footprint_m2": round(st["footprint_area_m2"], 1), "footprint": r(st["footprint_local_xy"], 2),
        "buffer_sweep": [[s["buffer_m"], s["points"]] for s in st["buffer_sweep"]],
        "xyz": r(xyz, 2), "labels": model.seg.labels.tolist(),
        "facets": [{k: f[k] for k in ("id", "points", "pitch_deg", "pitch_rise_per_12", "facing", "area_m2", "plan_area_m2",
                                      "confidence", "area_error_p90_pct", "outline", "area_not_counted_m2")}
                   | {"vertices": r(f["vertices_local_m"]), "holes": [r(x) for x in f.get("holes_local_m", [])]}
                   for f in rep["facets"]],
        "edges": [{"id": e["id"], "type": e["type"], "length_m": e["length_m"], "p0": r(e["p0_local_m"]),
                   "p1": r(e["p1_local_m"])} for e in rep["edges"]],
        "mesh": {"vertices": r(mesh["vertices"]), "faces": mesh["faces"].tolist(),
                 "facet_of_face": mesh["facet_of_face"].tolist()},
        "totals": rep["totals"], "quality": rep["quality"], "settings": rep["settings"],
    }


def excluded(h: str) -> dict:
    st = json.loads((ROOT / h / "stats.json").read_text())
    return {"id": h, "footprint_m2": round(st["footprint_area_m2"], 1),
            "buffer_sweep": [[s["buffer_m"], s["points"]] for s in st["buffer_sweep"]]}


def synthetic_area_errors(seeds: int = 20) -> list[dict]:
    """|area error| of every facet of the synthetic roofs through the full pipeline (report areas),
    matched to the true facet that most of its points belong to."""
    out = []
    for make in (gable, hip, cross_gable, dormer):
        for seed in range(seeds):
            roof = make(seed=seed)
            model = build_roof(roof.xyz, density=33.0, seed=seed)
            rep = measurement_report(model, {"house_id": make.__name__})
            facets = {f["id"]: f for f in rep["facets"]}
            labels = model.seg.labels
            for t, true_area in enumerate(roof.true_areas):
                hits = labels[(roof.labels == t) & (labels >= 0)]
                if len(hits) == 0:
                    continue
                s = int(np.bincount(hits).argmax())
                f = facets.get(s)
                # A facet RANSAC missed: its points sit in another facet's segment. Not a match.
                if f is None or np.bincount(roof.labels[labels == s]).argmax() != t:
                    continue
                out.append({"roof": make.__name__, "seed": seed, "true_m2": round(true_area, 2),
                            "error_pct": round(100 * (f["area_m2"] / true_area - 1), 2),
                            "fallback": f["outline"] != "plane and line intersections"})
    return out


def main() -> None:
    sweep = list(csv.DictReader((ROOT / "sweeps/stage2_sweep.csv").open()))
    joint = list(csv.DictReader((ROOT / "sweeps/stage2_joint_ransac.csv").open()))
    hand = Path("validation/hand_measurements.csv")
    data = {
        "houses": [house(h) for h in HOUSES],
        "excluded": [excluded("delft_04")],
        "validation": json.loads((ROOT / "validation.json").read_text()),
        "hand_measurements": len(list(csv.DictReader(hand.open()))) if hand.exists() else 0,
        "sweep": [{k: (v if k in ("setting", "method", "house", "reference_ok") else float(v)) for k, v in row.items()}
                  for row in sweep],
        "joint": [{k: float(v) for k, v in row.items()} for row in joint],
        "synthetic_area": synthetic_area_errors(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")
    errs = data["synthetic_area"]
    for lo, hi in ((0, 10), (10, 25), (25, np.inf)):
        e = np.abs([x["error_pct"] for x in errs if lo <= x["true_m2"] < hi])
        print(f"  synthetic facets {lo}-{hi} m²: n {len(e)}, |error| median {np.median(e):.1f}%, "
              f"p90 {np.percentile(e, 90):.1f}%, max {e.max():.1f}%")


if __name__ == "__main__":
    main()
