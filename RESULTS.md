# Results

Measured numbers only. Each row comes from a script output; the source is noted.

## Stage 1: load and isolate

Source: `scripts/stage1_isolate.py`. Density is per m² of plan area, measured on building
points inside the walls. "Non-building inside walls" is total / more than 0.5 m from the walls.

| house | footprint m² | buffer m | roof pts | nominal pts/m² | local pts/m² (r = 1 m) | median NN spacing m | non-building inside walls |
|---|---|---|---|---|---|---|---|
| delft_01 | 171.5 | 1.0 | 6148 | 33.5 | 33.1 | 0.089 | 539 / 76 |
| delft_02 | 92.3 | 1.0 | 4665 | 49.5 | 48.7 | 0.069 | 793 / 77 |
| delft_03 | 119.8 | 1.0 | 4163 | 31.6 | 30.6 | 0.093 | 562 / 92 |
| delft_05 | 85.3 | 1.0 | 3052 | 34.0 | 33.1 | 0.088 | 778 / 296 |

Across the 4 kept houses: nominal density 31.6–49.5 pts/m² (median 33.8), median NN spacing
0.069–0.093 m. Local and nominal density agree within 3.5% on every house.

Excluded:

| house | footprint m² | reason | sweep 0 / 1 / 2 / 3 m |
|---|---|---|---|
| delft_04 | 123.0 | sweep never plateaus: building-class roof continues past the outline (structure not in BAG) | 5261 / 5834 / 6030 / 6210 |

Houses come from `scripts/find_candidates.py`: of 700 BAG buildings on tile (withheld), 5 are
single-dwelling residential, 70–300 m², ≥ 3 m from any other building. All 5 were run.

All on AHN5 2023, GeoTiles sub-tile (withheld). Roof eave extent past the BAG outline (where
the sweep plateaus): 01: ≤ 0.75 m, 02: ≤ 0.5 m, 03: ≤ 0.75 m, 05: ≤ 0.5 m.

| house | BAG id | address | built | roof (from preview) |
|---|---|---|---|---|
| delft_01 | <bag-id> | (address withheld) | 1910s | hipped with flat top + flat annex |
| delft_02 | <bag-id> | (address withheld) | 1970s | two flat levels (6.4 m, 3.4 m): pitch-0 control |
| delft_03 | <bag-id> | (address withheld) | 1920s | gable + low flat extension + small lean-to |
| delft_04 | <bag-id> | (address withheld) | 1950s | excluded |
| delft_05 | <bag-id> | (address withheld) | 2010s | hipped |

Notes:

- delft_01 buffer sweep (building points): 0 m: 5749 · 0.25 m: 6142 · 0.5 m: 6147 · 0.75 m: 6148 · 1.0–3.0 m: 6148.
  Points stop by 0.5 m, so the eave overhang is ≤ 0.5 m; 1.0 m sits on the plateau, and
  neither shed (3.7 m, 7.2 m away) is building-class.
- delft_01: 6.5% of the house's building points (399 / 6148) lie outside the walls, in the eaves.
- delft_01: 463 of the 539 non-building points inside the walls are within 0.5 m of the wall line
  (facade and ground at the foot of the wall). Of the 76 interior ones, 19 are above the
  ridge: an unclassified chimney, which also leaves a hole in the building-class roof.
- delft_01: 17 building-class points are below z = 1.5 m along the north wall (ground or wall labelled
  building). Kept for now; stage 2 has to reject them.
- delft_05: 296 non-building points 0.5–1.5 m inside the outline on two sides, z 0–3 m,
  below the roof (≥ 4.2 m): facade or ground seen under the roof edge. They don't touch the roof.
- delft_02, 03: 12 and 18 non-building points above the roof's highest building point
  (chimneys or antennas).
- No dormer or L-plan-with-valley house on this tile; that needs another tile.

## Stage 2: plane segmentation

### Synthetic checks: plane fit on known facets

Source: `tests/synthetic_roofs.py` + `roofrecon.planes`, 20 seeds, 33 pts/m², 3 cm vertical noise,
facet membership given (no segmentation yet). Area = convex hull in the plane.

| roof | facets | pitch error ° (median / max abs) | area error % (median / p5 / p95) |
|---|---|---|---|
| gable 35° | 40 | +0.013 / 0.078 | −1.10 / −2.09 / −0.48 |
| hip 30° | 80 | +0.010 / 0.151 | −1.52 / −3.19 / −0.64 |

Area always reads low: the outermost points sit inside the true edge. The bias is larger on
triangular hip ends (more edge per m²).

### Synthetic checks: sequential RANSAC

Source: `roofrecon.segment`, 10 seeds per row, 33 pts/m², 3 cm vertical noise, min_points 30,
connect_radius 0.3 m, 500 iterations. "raw" = `sequential_ransac`; "clean" = `segment_ransac`
(reassign → merge coplanar → reassign). Mislabelled = points not in their segment's majority facet.

| roof | threshold m | right plane count raw → clean | mislabelled % median raw → clean (clean max) | coverage % (clean, median) |
|---|---|---|---|---|
| gable 35° | 0.05 | 100 → 100 | 0.43 → 0.36 (0.53) | 95.2 |
| gable 35° | 0.10 | 100 → 100 | 0.99 → 0.34 (0.53) | 99.5 |
| gable 35° | 0.20 | 100 → 100 | 2.43 → 0.34 (0.53) | 99.6 |
| gable 10° | 0.05 | 100 → 100 | 1.84 → 1.61 (1.92) | 90.8 |
| gable 10° | 0.10 | 100 → 100 | 5.00 → 2.01 (2.33) | 99.6 |
| gable 10° | 0.20 | 100 → 100 | 14.15 → 3.13 (3.56) | 99.7 |
| hip 30° | 0.05 | 90 → 100 | 1.89 → 1.49 (1.92) | 93.5 |
| hip 30° | 0.10 | 80 → 100 | 3.58 → 1.54 (1.95) | 99.4 |
| hip 30° | 0.20 | 80 → 100 | 7.86 → 1.58 (1.95) | 99.2 |
| hip 15° | 0.05 | 90 → 90 | 3.61 → 3.25 (3.94) | 91.5 |
| hip 15° | 0.10 | 90 → 100 | 7.32 → 3.56 (4.25) | 99.4 |
| hip 15° | 0.20 | 30 → 30 | 57.81 → 56.31 (62.99) | 91.1 |

- 5 cm threshold: coverage 91–95%; perpendicular noise is 3 cm × cos(pitch), so a few % of
  points fall outside ±5 cm.
- Ridge stealing grows with threshold and shrinks with pitch (strip width ≈ threshold / sin(2·pitch)).
- Hip 15° at 0.20 m fails: a near-horizontal slab through the ridge area outscores every real
  facet. Clean-up can't fix a wrong first plane; needs per-point normals.
- Raw hip errors at 0.10–0.20 m were extra planes: hip-end tips cut off by strips stolen
  along both hips. Reassigning before merging reconnects them.
