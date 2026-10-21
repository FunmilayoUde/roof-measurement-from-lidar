"""Stage 4: write the measurement report for one house.

    uv run python scripts/stage4_report.py --house-id delft_03

Reads data/interim/<house-id>/{roof.ply, stats.json, footprint.geojson}, builds the roof
model with the chosen stage 2/3 settings, and writes data/interim/<house-id>/report.json.
Coordinates in the report are local (stage 1 recentering); add source.offset_xyz to get
EPSG:28992 / NAP.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import open3d as o3d

from roofrecon.measure import measurement_report
from roofrecon.pipeline import build_roof


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    house = args.root / args.house_id
    xyz = np.asarray(o3d.io.read_point_cloud(str(house / "roof.ply")).points)
    stats = json.loads((house / "stats.json").read_text())
    footprint = json.loads((house / "footprint.geojson").read_text())
    meta = {
        "house_id": args.house_id,
        "bag_id": footprint.get("properties", {}).get("identificatie"),
        "tile": Path(stats["tile"]["path"]).name,
        "crs": stats["tile_crs"],
        "offset_xyz": stats["offset_xyz"],
        "footprint_area_m2": round(stats["footprint_area_m2"], 1),
    }
    rep = measurement_report(build_roof(xyz, stats["density_inside_walls"]["nominal_pts_per_m2"], seed=args.seed), meta)
    (house / "report.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False))

    t = rep["totals"]
    print(f"{args.house_id}: {t['facets']} facets, {t['roof_area_m2']} m² = {t['roof_area_ft2']:.0f} ft² = {t['squares']} squares")
    print(f"  {'facet':>5} {'pitch':>12} {'facing':>6} {'area m²':>8}  outline")
    for f in sorted(rep["facets"], key=lambda f: -f["area_m2"]):
        print(f"  {f['id']:>5} {f['pitch_deg']:>5}° ({f['pitch_rise_per_12']:>4}/12) {f['facing'] or 'flat':>6} {f['area_m2']:>8}  {f['outline']}")
    print("  lengths: " + ", ".join(f"{k} {v} m ({t['length_ft_by_type'][k]} ft)" for k, v in t["length_m_by_type"].items()))
    for note in rep["quality"]["notes"]:
        print(f"  note: {note}")


if __name__ == "__main__":
    main()
