# Roof Measurement from LiDAR Alone

Python package: `roofrecon`.

Roof reconstruction from open aerial LiDAR: facet polygons with pitch and true area, a typed
wireframe (ridge, hip, valley, eave, rake), and a measurement report. Numbers go in
[RESULTS.md](RESULTS.md).

## Setup

```sh
uv sync                 # Python 3.12 (Open3D has no wheels for newer versions yet)
uv run pytest
```

## Running the pipeline on one house

```sh
uv run python scripts/stage1_isolate.py --laz data/raw/<tile>.LAZ --bag-id <BAG pand id> --house-id <name> --crs EPSG:28992
uv run python scripts/stage2_segment.py --house-id <name> --method ransac   # planes, previews, normal map
uv run python scripts/stage3_topology.py --house-id <name>                 # facet polygons, typed edges
uv run python scripts/stage4_report.py --house-id <name>                   # report.json
uv run python scripts/export_model.py --house-id <name>                   # model.obj / .ply, edges, preview
uv run python scripts/view_model.py --house-id <name>                     # interactive 3D view over the points
uv run python scripts/stage5_validate.py                                   # footprint checks, hand measurements
uv run python scripts/make_report_data.py                                 # data for the write-up, docs/report/
```

Finding houses: `scripts/scout_tiles.py` scores AHN sub-tiles around a place before downloading;
`scripts/find_candidates.py` lists detached houses in a tile or bbox. Settings were chosen with
`scripts/stage2_sweep.py`. Example output: `docs/examples/`. Project write-up: `docs/report/index.html` (serve the folder
locally, e.g. `python -m http.server -d docs/report`; it loads `report_data.json`).

## Layout

```
src/roofrecon/
  las_io.py      chunked LAZ reading with a bbox filter; header summary
  footprint.py   OSM building lookup (Overpass) and reprojection
  isolate.py     class filter, footprint crop, buffer sweep, recentering, density
  planes.py      plane fit, pitch, azimuth, areas, alpha shape
  normals.py     per-point normals by local PCA
  segment.py     sequential RANSAC, region growing, clean-up
  topology.py    plane intersections, corners, ridge/hip/valley/break
  polygons.py    facet polygons with eave/rake/top/perimeter sides
  pipeline.py    stages 2-3 with the chosen settings
  measure.py     the stage 4 report
  export.py      the 3D roof model (mesh + edges)
  validate.py    stage 5 checks
  report.py      per-plane tables for stage 2
  viz.py         previews, normal map, .ply writing
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
- **Buffer beyond the footprint.** On every house so far the roof extends past the BAG
  outline, by 0.25–0.75 m (RESULTS.md), so treat the outline as roughly the walls. Crop at the walls and
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
- **Density, measured in plan view, on building points inside the walls.** The roof
  covers that whole region, so the area isn't diluted by the empty part of the buffer
  (dividing by the buffered area gave 25.1 instead of 33.5 on delft_01). `nominal` =
  points / footprint area. `local` = median neighbour count within 1 m / π. Nearest-neighbour spacing is reported on
  its own because converting it to density needs an assumption about the sampling
  pattern.
