"""Stage 1: isolate one house's roof points from a LAZ tile.

    uv run python scripts/stage1_isolate.py --laz data/raw/<tile>.LAZ \
        --lat 52.0 --lon 4.36 --house-id delft_01

Writes data/interim/<house-id>/{footprint.geojson, roof.ply, stats.json, preview.png}
and prints the row to paste into RESULTS.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from roofrecon.footprint import fetch_osm_building, load_footprint
from roofrecon.isolate import (
    ASPRS_CLASSES,
    buffer_sweep,
    classification_histogram,
    crop_mask,
    density_stats,
    recenter,
)
from roofrecon.las_io import header_summary, read_points
from roofrecon.viz import save_snapshot, write_ply

SWEEP = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--laz", required=True, type=Path)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--footprint", type=Path, help="GeoJSON Feature; fetched from OSM if absent and --lat/--lon given")
    ap.add_argument("--footprint-crs", default="EPSG:4326")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    ap.add_argument("--crs", help="Tile CRS if the LAZ header lacks one, e.g. EPSG:28992 (AHN) or a UTM zone")
    ap.add_argument("--classes", default="6", help="Comma-separated ASPRS classes to keep (6 = building)")
    ap.add_argument("--buffer", type=float, default=1.0, help="Metres beyond the wall footprint; see the sweep")
    ap.add_argument("--out-root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    out = args.out_root / args.house_id
    out.mkdir(parents=True, exist_ok=True)

    header = header_summary(args.laz)
    tile_crs = args.crs or header["crs"]
    if tile_crs is None:
        raise SystemExit("LAZ header has no CRS; pass --crs (AHN is EPSG:28992)")

    fp_path = args.footprint or out / "footprint.geojson"
    if not fp_path.exists():
        if args.lat is None or args.lon is None:
            raise SystemExit("Need --footprint or --lat/--lon")
        fetch_osm_building(args.lat, args.lon, fp_path)
    footprint = load_footprint(fp_path, src_crs=args.footprint_crs, dst_crs=tile_crs)

    # Read a generous box so the sweep can see past the eaves.
    xmin, ymin, xmax, ymax = footprint.buffer(max(SWEEP) + 1.0).bounds
    pts = read_points(args.laz, bbox=(xmin, ymin, xmax, ymax))
    if len(pts) == 0:
        raise SystemExit("No points in the footprint's bbox: wrong tile or wrong CRS?")
    xy = pts.xyz[:, :2]

    # What the classifier says lives inside the walls. Non-building points here are either
    # genuine (trees overhanging the roof) or misclassified roof, which later caps coverage.
    inside_walls = crop_mask(xy, footprint, 0.0)
    hist_inside = classification_histogram(pts.classification[inside_walls])

    keep_classes = [int(c) for c in args.classes.split(",")]
    building = pts.subset(np.isin(pts.classification, keep_classes))
    sweep = buffer_sweep(building.xyz[:, :2], footprint, SWEEP)

    roof = building.subset(crop_mask(building.xyz[:, :2], footprint, args.buffer))
    local, offset = recenter(roof.xyz)

    crop_region = footprint.buffer(args.buffer, join_style="mitre")
    stats = {
        "house_id": args.house_id,
        "tile": header,
        "tile_crs": tile_crs,
        "offset_xyz": offset.tolist(),
        "footprint_area_m2": footprint.area,
        "footprint_local_xy": (np.asarray(footprint.exterior.coords) - offset[:2]).tolist(),
        "classes_kept": keep_classes,
        "buffer_m": args.buffer,
        "class_histogram_inside_walls": {f"{k} {ASPRS_CLASSES.get(k, '?')}": v for k, v in hist_inside.items()},
        "buffer_sweep": sweep,
        "roof_points": len(roof),
        "z_range_m": [float(local[:, 2].min()), float(local[:, 2].max())],
        "density_crop_region": density_stats(local[:, :2], crop_region.area),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    write_ply(local, out / "roof.ply")
    save_snapshot(
        local,
        out / "preview.png",
        footprint_xy=np.asarray(footprint.exterior.coords) - offset[:2],
        title=f"{args.house_id}: {len(roof)} pts, buffer {args.buffer} m",
    )

    d = stats["density_crop_region"]
    frac_building = hist_inside.get(6, 0) / max(sum(hist_inside.values()), 1)
    print("\nclass histogram inside walls:", stats["class_histogram_inside_walls"])
    print("buffer sweep (building points):")
    for r in sweep:
        print(f"  {r['buffer_m']:>5.2f} m  {r['points']:>7d}  (+{r['added'] if r['added'] is not None else '-'})")
    print("\nRESULTS.md row:")
    print(
        f"| {args.house_id} | {footprint.area:.1f} | {args.buffer} | {len(roof)} | "
        f"{d['nominal_pts_per_m2']:.1f} | {d['local_pts_per_m2']:.1f} | {d['median_nn_spacing_m']:.3f} | "
        f"{frac_building:.3f} |"
    )


if __name__ == "__main__":
    main()
