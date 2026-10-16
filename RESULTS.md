# Results

Measured numbers only. Each row comes from a script output; the source is noted.

## Stage 1: load and isolate

Source: `scripts/stage1_isolate.py`. Density is per m² of plan area, measured on building
points inside the walls. "Non-building inside walls" is total / more than 0.5 m from the walls.

| house | footprint m² | buffer m | roof pts | nominal pts/m² | local pts/m² (r = 1 m) | median NN spacing m | non-building inside walls |
|---|---|---|---|---|---|---|---|
| delft_01 | 171.5 | 1.0 | 6148 | 33.5 | 33.1 | 0.089 | 539 / 76 |

delft_01: BAG <bag-id> (Delft, built (decade withheld)), AHN5 2023, GeoTiles sub-tile (withheld).

- Buffer sweep (building points): 0 m: 5749 · 0.25 m: 6142 · 0.5 m: 6147 · 0.75 m: 6148 · 1.0–3.0 m: 6148.
  Points stop by 0.5 m, so the eave overhang is ≤ 0.5 m; 1.0 m sits on the plateau, and
  neither shed (3.7 m, 7.2 m away) is building-class.
- 6.5% of the house's building points (399 / 6148) lie outside the walls, in the eaves.
- 463 of the 539 non-building points inside the walls are within 0.5 m of the wall line
  (facade and ground at the foot of the wall). Of the 76 interior ones, 19 are above the
  ridge: an unclassified chimney, which also leaves a hole in the building-class roof.
- 17 building-class points are below z = 1.5 m along the north wall (ground or wall labelled
  building). Kept for now; stage 2 has to reject them.
