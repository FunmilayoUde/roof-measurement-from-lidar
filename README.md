# roofrecon

Roof reconstruction from open aerial LiDAR: facet polygons with pitch and true area, a typed
wireframe (ridge, hip, valley, eave, rake), and a measurement report. Numbers go in
[RESULTS.md](RESULTS.md).

## Setup

```sh
uv sync                 # Python 3.12 (Open3D has no wheels for newer versions yet)
uv run pytest
```

## Layout

```
src/roofrecon/
  las_io.py      chunked LAZ reading with a bbox filter; header summary
  footprint.py   OSM building lookup (Overpass) and reprojection
  isolate.py     class filter, footprint crop, buffer sweep, recentering, density
  viz.py         matplotlib snapshot, Open3D conversion, .ply writing
scripts/
  stage1_isolate.py   tile + house -> data/interim/<id>/{roof.ply, stats.json, preview.png}
  view.py             Open3D viewer
tests/            synthetic scenes with known answers
data/             not in git: raw/ for tiles, interim/ for per-house outputs
```

## Stage 1: load and isolate

### Getting a tile (AHN, Netherlands)

1. Pick a **detached** house (see "Choosing houses" below) and note its lat/lon.
2. Find the AHN tile that covers it. GeoTiles (TU Delft) splits AHN into 1 x 1.25 km
   sub-tiles, a few hundred MB each instead of the 4+ GB national sheets. Tile pages
   are at `https://geotiles.citg.tudelft.nl/tiles/html/<TILE>.html`. Download the
   AHN4 (or AHN5 where available) `.LAZ` into `data/raw/`.
3. Run:

```sh
uv run python scripts/stage1_isolate.py --laz data/raw/<tile>.LAZ \
    --lat <lat> --lon <lon> --house-id <name> --crs EPSG:28992
uv run python scripts/view.py data/interim/<name>/roof.ply
```

Open `roof.ply` in CloudCompare as well. Because the coordinates are recentered, it
won't ask for a global shift.

### Choosing houses

- **Detached, not terraced.** Most Dutch single-family homes are row houses whose roofs
  run continuously across party walls. The footprint splits the building, but no
  geometric edge does. Stages 2 and 3 would then produce facets that end at a wall
  you can't see in the points. Begin with detached houses. Terraced ones are a later,
  separate problem.
- Vary the roof types: gable, hip, and at least one with a dormer or an L-shaped plan.
- Avoid heavy tree cover over the roof. Points under the canopy are missing, not
  misplaced.

### Why each step exists

- **Class filter (6 = building).** The provider's classifier already separates ground,
  vegetation, and buildings. We trust it for now and measure how much it gets wrong:
  the class histogram inside the walls shows what *isn't* labelled building, and
  anything there caps stage-2 coverage.
- **Buffer beyond the footprint.** OSM and BAG footprints trace the *walls*, and the
  roof overhangs them by the eave depth (often 0.3 to 0.6 m). Crop at the walls and
  you cut the eaves off, which biases area low and moves the eave edges inward. The
  buffer sweep counts building points as the crop grows. The count should rise up to
  the eave, stay flat, then rise again at the next building. Pick a buffer on the
  plateau.
- **Mitred buffer.** A round buffer cuts off square eave corners. Mitred joins keep
  them.
- **Recentering.** Projected coordinates are large. float32 (used by most PLY
  consumers and GPUs) resolves about 3 cm at y ≈ 450 000 m, and about 50 cm at a UTM
  northing of 4.5 million m. Subtracting a whole-metre offset brings that to
  sub-micron. The offset goes in `stats.json` so results can be mapped back.
- **Density, measured in plan view.** `nominal` = points / crop area. `local` =
  median neighbour count within 1 m / π. Nearest-neighbour spacing is reported on
  its own because converting it to density needs an assumption about the sampling
  pattern.
