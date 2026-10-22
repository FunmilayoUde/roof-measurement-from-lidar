"""Stage 1: isolate one house's roof points from a LAZ tile.

    uv run python scripts/stage1_isolate.py --laz data/raw/<tile>.LAZ \
        --lat 52.0 --lon 4.36 --house-id delft_01      # OSM footprint
    uv run python scripts/stage1_isolate.py --laz data/raw/<tile>.LAZ \
        --bag-id <BAG pand id> --house-id delft_01     # Dutch BAG footprint

Writes data/interim/<house-id>/{footprint.geojson, roof.ply, stats.json, preview.png}
and prints the row to paste into RESULTS.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from roofrecon.footprint import fetch_bag_building, fetch_osm_building, load_footprint
from roofrecon.isolate import (
    ASPRS_CLASSES,
    buffer_sweep,
    classification_histogram,
    crop_mask,
    density_stats,
    non_building_inside_walls,
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
    ap.add_argument("--bag-id", help="Dutch BAG pand id; footprint fetched from PDOK in EPSG:28992")
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
    if args.bag_id:
        args.footprint_crs = "EPSG:28992"
    if not fp_path.exists():
        if args.bag_id:
            fetch_bag_building(args.bag_id, fp_path)
        elif args.lat is not None and args.lon is not None:
            fetch_osm_building(args.lat, args.lon, fp_path)
        else:
            raise SystemExit("Need --footprint, --bag-id, or --lat/--lon")
    footprint = load_footprint(fp_path, src_crs=args.footprint_crs, dst_crs=tile_crs)

    # Read a generous box so the sweep can see past the eaves.
    xmin, ymin, xmax, ymax = footprint.buffer(max(SWEEP) + 1.0).bounds
    pts = read_points(args.laz, bbox=(xmin, ymin, xmax, ymax))
    if len(pts) == 0:
        raise SystemExit("No points in the footprint's bbox: wrong tile or wrong CRS?")
    xy = pts.xyz[:, :2]

    # What the classifier says lives inside the walls. Non-building points at the wall line
    # are facade/ground hits; interior ones are unlabelled roof objects or trees.
    inside_walls = crop_mask(xy, footprint, 0.0)
    hist_inside = classification_histogram(pts.classification[inside_walls])
    nonb = non_building_inside_walls(xy, pts.classification, footprint)

    keep_classes = [int(c) for c in args.classes.split(",")]
    building = pts.subset(np.isin(pts.classification, keep_classes))
    sweep = buffer_sweep(building.xyz[:, :2], footprint, SWEEP)

    roof = building.subset(crop_mask(building.xyz[:, :2], footprint, args.buffer))
    local, offset = recenter(roof.xyz)

    # Density inside the walls: a region the roof fully covers, so the area is not diluted by
    # empty buffer. Eaves outside the walls are left out of the estimate, not the output.
    in_walls_local = local[crop_mask(roof.xyz[:, :2], footprint, 0.0)]
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
        "non_building_inside_walls": {k: v for k, v in nonb.items() if k != "mask_interior"},
        "non_building_interior_xyz_local": (pts.xyz[nonb["mask_interior"]] - offset).round(3).tolist(),
        "z_range_m": [float(local[:, 2].min()), float(local[:, 2].max())],
        "density_inside_walls": density_stats(in_walls_local[:, :2], footprint.area),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    write_ply(local, out / "roof.ply")
    save_snapshot(
        local,
        out / "preview.png",
        footprint_xy=np.asarray(footprint.exterior.coords) - offset[:2],
        title=f"{args.house_id}: {len(roof)} pts, buffer {args.buffer} m",
    )

    d = stats["density_inside_walls"]
    print("\nclass histogram inside walls:", stats["class_histogram_inside_walls"])
    print(
        f"non-building inside walls: {nonb['total']} total, "
        f"{nonb['interior']} more than {nonb['wall_band_m']} m from the walls"
    )
    print("buffer sweep (building points):")
    for r in sweep:
        print(f"  {r['buffer_m']:>5.2f} m  {r['points']:>7d}  (+{r['added'] if r['added'] is not None else '-'})")
    print("\nRESULTS.md row:")
    print(
        f"| {args.house_id} | {footprint.area:.1f} | {args.buffer} | {len(roof)} | "
        f"{d['nominal_pts_per_m2']:.1f} | {d['local_pts_per_m2']:.1f} | {d['median_nn_spacing_m']:.3f} | "
        f"{nonb['total']} / {nonb['interior']} |"
    )


if __name__ == "__main__":
    main()
