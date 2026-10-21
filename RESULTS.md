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

### Synthetic checks: normals

Source: `roofrecon.normals.estimate_normals`, gable, 33 pts/m², 3 cm noise. Median normal error (°)
by plan distance from the ridge.

| radius m | neighbours | 0–0.1 m | 0.2–0.3 m | 0.4–0.6 m | > 1.5 m |
|---|---|---|---|---|---|
| 0.3 | 8 | 26.5 | 5.6 | 5.2 | 5.32 |
| 0.5 | 21 | 27.9 | 8.1 | 1.6 | 1.49 |
| 0.75 | 47 | 31.5 | 15.6 | 2.6 | 0.60 |
| 1.0 | 83 | 32.3 | 20.6 | 8.7 | 0.34 |

(35° gable; 15° gable follows the same pattern with roughly half the error at the ridge.)
Larger radius: less noise in the interior, wider blurred band at the ridge.

### Synthetic checks: RANSAC vs region growing

Same setup as the RANSAC table. Region growing: normal radius 0.5 m, angle 8°, grow radius 0.3 m,
seeds only where surface variation ≤ 2 × median; then merge → absorb → reassign → merge → reassign.

| roof | threshold m | right count RANSAC / RG | mislabelled % median RANSAC / RG | coverage % RANSAC / RG |
|---|---|---|---|---|
| gable 35° | 0.05 | 100 / 100 | 0.36 / 0.32 | 95.2 / 95.9 |
| gable 35° | 0.10 | 100 / 100 | 0.34 / 0.34 | 99.5 / 99.5 |
| gable 35° | 0.20 | 100 / 100 | 0.34 / 0.34 | 99.6 / 99.5 |
| gable 10° | 0.05 | 100 / 100 | 1.61 / 1.54 | 90.8 / 92.2 |
| gable 10° | 0.10 | 100 / 100 | 2.01 / 1.83 | 99.6 / 99.8 |
| gable 10° | 0.20 | 100 / 100 | 3.13 / 1.84 | 99.7 / 99.9 |
| hip 30° | 0.05 | 100 / 80 | 1.49 / 1.49 | 93.5 / 94.7 |
| hip 30° | 0.10 | 100 / 100 | 1.54 / 1.51 | 99.4 / 99.4 |
| hip 30° | 0.20 | 100 / 100 | 1.58 / 1.52 | 99.2 / 99.6 |
| hip 15° | 0.05 | 90 / 100 | 3.25 / 3.16 | 91.5 / 92.8 |
| hip 15° | 0.10 | 100 / 90 | 3.56 / 3.39 | 99.4 / 99.7 |
| hip 15° | 0.20 | 30 / 90 | 56.31 / 3.42 | 91.1 / 99.9 |

Region growing's 4 wrong counts (of 120 runs):
- hip 30° at 0.05 m, two runs: one extra 36-point plane of blended normals along a hip line
  (facing 231°, between south and west); one hip end split in two.
- hip 15° at 0.10 and 0.20 m, the same seed: one extra 51-point plane at 2.5° on top of the
  ridge, where the blend of two opposite 15° slopes is nearly flat.
The extra planes are ≈ 1–1.5 m², just above min_points (30); real houses should set that.

### Real houses: expected planes (written before running stage 2)

Read off the stage-1 previews (plan + two elevations) only. A range where the preview is ambiguous.

| house | expected planes | reasoning |
|---|---|---|
| delft_01 | 6 (6–7) | main block: hip with a flat top = 4 slopes + 1 flat; flat annex at 2.8 m = 1. Unsure about the link between annex and main block. |
| delft_02 | 2 | two flat levels, 6.4 m and 3.4 m. Pitch should be ~0° on both. |
| delft_03 | 4 (3–5) | gable = 2; low flat extension at ~3.7 m = 1; small lean-to on one end = 1, but it is small and may fall under min_points. |
| delft_05 | 4 (4–6) | hip; the plan has notched corners and two high spots, so possibly a split ridge with extra facets. |

Not expected as planes: chimneys (mostly unlabelled), the 17 low building-class points on delft_01.
Any plane steeper than ~75° is a wall, not a roof facet.

Pipeline change after the synthetic tables above: both methods now end with "keep each plane's
largest connected piece → drop planes under min_points → absorb by adjacency". Rerun on the
synthetic set: plane counts unchanged in every row; mislabelled within ±0.5 points and
coverage within ±0.8 points of the tables above.

### Real houses: first run

Source: `scripts/stage2_segment.py`, threshold 0.10 m, min_area 1 m² (× house density), connect /
grow radius 0.3 m; region growing: normal radius 0.5 m, angle 8°. Area = convex hull in the plane.
Walls (> 75°): none found.

| house | expected | RANSAC planes | RG planes | coverage R / RG | median rms m R / RG | roof true area m² R / RG |
|---|---|---|---|---|---|---|
| delft_01 | 6 (6–7) | 9 | 15 | 0.828 / 0.899 | 0.029 / 0.025 | 224.3 / 240.4 |
| delft_02 | 2 | 2 | 3 | 0.952 / 0.944 | 0.032 / 0.025 | 93.3 / 94.1 |
| delft_03 | 4 (3–5) | 5 | 5 | 0.964 / 0.963 | 0.027 / 0.027 | 162.9 / 163.2 |
| delft_05 | 4 (4–6) | 4 | 7 | 0.894 / 0.916 | 0.035 / 0.031 | 100.7 / 102.7 |

- delft_01, prediction missed: it's a mansard (four sides at 55–56°, both methods), and the
  "flat" top isn't planar. A cross-section shows two smooth humps ~0.25 m high and ~2 m wide
  along the roof, plus a smaller bump and the chimney. Both methods slice the top into stripes.
- delft_02: RANSAC matches (0.4° and 0.2°). RG adds a 122-point, 2.8 m², 12.6° plane, at a small
  round feature visible in the normal map near the roof edge.
- delft_03: within range, and the two methods agree. Matches the 5 colour patches in the normal
  map. Gable halves 55.4 / 56.0 m² at 42.7° / 44.1° (RANSAC).
- delft_05: count matches for RANSAC, but unresolved. The normal map shows the SE slope's facing
  changing gradually, plus small corner patches. RANSAC folds the corners into the SE slope; RG
  separates them but also splits the NW slope (301° / 314°).
- Before stray removal, reassignment by distance alone handed far-away coplanar points to planes:
  delft_03's south gable half read 84.3 m² instead of 55.9 m².

### Parameter sweeps on real houses

Source: `scripts/stage2_sweep.py` (172 one-at-a-time runs, RANSAC seed 0; then a RANSAC
threshold × min_area grid over seeds 0–9). Figures: `docs/figures/sweep_*.png`.
References: delft_02 must give exactly 2 flat planes; delft_03 exactly 5, gable halves within 5%
in area. delft_01 and delft_05 are reported but don't decide settings.

One-at-a-time (other settings at threshold 0.10, min_area 1.0, normal radius 0.5, angle 8°):

| setting | RANSAC passes both references | region growing passes both references |
|---|---|---|
| threshold 0.04–0.25 m | 0.08–0.25 | never (delft_02 always gets a 3rd plane) |
| min_area 0.25–4 m² | 1.0–2.0 | never (delft_03 needs exactly 1.0, delft_02 needs ≥ 2.0) |
| normal radius 0.3–1.0 m | n/a | never (delft_03 passes only at 0.5) |
| angle 4–20° | n/a | never (delft_03 passes at 8–20°) |

RANSAC joint grid, fraction of 10 seeds passing both references:

| threshold \ min_area | 0.75 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|
| 0.08 | 0% | 30% | 100% | 90% |
| 0.10 | 100% | 100% | 100% | 90% |
| 0.12 | 100% | 100% | 100% | 10% |
| 0.15 | 100% | 100% | 100% | 0% |
| 0.20 | 100% | 100% | 30% | 0% |

The seed-0 sweep made min_area 2.0 look safe; over 10 seeds it mostly isn't.

**Chosen: RANSAC, threshold 0.12 m, min_area 1.0 m²**, the centre of the all-seeds-pass region.
Seed stability there (10 seeds):

| house | planes | coverage | roof true area m² |
|---|---|---|---|
| delft_02 | 2 (10/10) | 0.971–0.973 | 94.3–94.5 |
| delft_03 | 5 (10/10) | 0.973–0.974 | 163.9–164.0 |
| delft_05 | 4 (10/10) | 0.931–0.940 | 101.5–102.2 |
| delft_01 | 6 (2), 8 (5), 9 (2), 11 (1) | 0.844–0.916 | 216.5–241.0 |

Trends:
- Coverage rises with threshold; on the references it is ≥ 0.97 from 0.12 m.
- Median plane rms rises roughly linearly with threshold (wider inlier band). On delft_05 RANSAC's
  rms climbs from 0.035 to 0.060 m between 0.10 and 0.25 m: separate slopes merging.
- Roof area on the references is flat from 0.08 m (delft_03 163–165 m², delft_02 93–96 m²).
- Below min_area 1.0, delft_03 gains a fragment; at ≥ 3.0 it loses its smallest real facet
  (2.5 m² in plan).
- Region growing: radius 0.3–0.4 m over-splits flat roofs (delft_02: 7 and 5 planes), ≥ 0.6 m
  loses delft_03's small facet; angle 4° over-splits (delft_02: 8 planes).
- Region growing's persistent 3rd plane on delft_02 is a ~2.8 m², 12.6° feature at a round object
  visible in the normal map. It is real, but it isn't a roof facet.

Decision: RANSAC carries into stage 3. Region growing has no setting that passes both references
on real data, though it fixed the shallow-hip slab failure on synthetic roofs; that case hasn't
appeared on a real house yet. The normal map stays as an independent check.

Final stage 2 run (RANSAC, threshold 0.12, min_area 1.0):

| house | planes | coverage | median rms m | roof true area m² |
|---|---|---|---|---|
| delft_02 | 2 | 0.973 | 0.037 | 94.5 |
| delft_03 | 5 | 0.973 | 0.029 | 164.0 |
| delft_05 | 4 | 0.933 | 0.037 | 101.5 |
| delft_01 | 8 | 0.844 | 0.032 | 223.3 |

## Stage 3: topology

### Synthetic checks: internal edges

Source: `roofrecon.topology.internal_edges` on RANSAC segmentations (threshold 0.12 m, 33 points),
20 seeds per roof, 33 pts/m², 3 cm noise. Endpoint error = distance to the nearer true endpoint.
"Snapped" endpoints are three-plane corners; the rest end where the contact points stop.

| roof | edges found / true | extra | type correct | snapped endpoint error cm (median / p90 / max) | other endpoint error cm (median / p90 / max) |
|---|---|---|---|---|---|
| gable 35° | 20 / 20 | 0 | 20 / 20 | none | 33.5 / 79.2 / 158.7 |
| hip 30° | 100 / 100 | 0 | 100 / 100 | 0.9 / 1.6 / 1.9 (n = 117) | 60.0 / 155.4 / 298.3 |
| cross-gable 35° | 80 / 80 | 0 | 80 / 80 | 0.5 / 1.0 / 1.4 (n = 44) | 21.8 / 113.9 / 279.8 |

- Three-plane corners are exact to ~1 cm. Endpoints without a third facet (at eave corners and
  gable ends) are 0.2–3 m short; step 2 replaces them with intersections with the outer edges.
- Two bugs found by the cross-gable (a T: a wing joining a main gable in two valleys):
  - Deciding which side of an edge a facet lies on from the facet's centroid fails for a notched
    facet, which has points on both sides of the edge line: 4 of 20 valleys read as breaks.
    Right at the line, labels from both facets mix (1 more seed). The side now comes from the
    facet's points alongside the edge, 0.2–1.5 m back from it.
  - When contact points stop short, both ends of an edge could snap to the same corner and
    collapse it (1 seed). Each corner now fixes at most one end.

### Synthetic checks: facet polygons

Source: `roofrecon.polygons.roof_polygons` on RANSAC segmentations (threshold 0.12 m, 33 points),
20 seeds per roof. Outline: alpha shape, radius 0.75 m. Area = polygon area in the plane.
Length error = total length of that edge type minus the true total.

| roof | facet area error % (median / min / max) | eave m | rake m | ridge m | hip m | valley m | fallbacks | "other" sides |
|---|---|---|---|---|---|---|---|---|
| gable 35° | −0.43 / −1.24 / +0.10 | −0.09 of 20.0 | −0.00 of 19.5 | −0.04 of 10.0 | | | 0 | 2 |
| hip 30° | −0.06 / −1.53 / +1.21 | −0.01 of 40.0 | | +0.00 of 4.0 | −0.01 of 24.4 | | 0 | 0 |
| cross-gable 35° | −0.54 / −2.43 / +0.47 | −0.11 of 32.0 | −0.05 of 24.4 | −0.07 of 20.0 | | −0.03 of 6.3 | 0 | 0 |

(length errors are medians over the 20 seeds)

Compared with stage 2's convex-hull area (median −1.1% gable, −1.5% hip, always low), the
polygon area is closer and no longer one-sided on the hip.

How it got there (each step measured on the same 20-seed set, median facet area error):
- shapely's concave_hull: its setting is relative to the longest edges, so no value works on all
  three roofs (−10% to −60% when tight, +7.7% on the notched facet when loose). Replaced by an
  alpha shape with an absolute radius: one piece, no holes for radius ≥ 0.35 m.
- Corners from line intersections, outer lines through the outline: −2.7% to −3.5%, with
  fallbacks. The outline rounds every corner, which leaves short diagonal sides.
- Dropping short and diagonal artefact sides, and merging consecutive parallel lines (they can't
  form a corner): fallbacks 0–4.
- Outer lines moved to the 2nd-outermost point alongside them (boundary points sit inside the true
  edge): −0.2% to −0.9%.
- Labelling the outline by which neighbouring facet it touches, instead of matching outline
  segments to edge lines by distance and angle: valleys from −1.6 m to −0.03 m.
- Stage 2 fix: keep a plane's minor pieces when within 1 m of its main piece. Dropping them left a
  97-point hole at a hip-end apex.

Remaining: 2 of 40 gable facets get an "other" side where missing points dent a corner (one rake
4.9 m short). The snap angle (10°) and artefact length (2 m) could be swept on real houses.

Later changes, rerun on the same 20-seed set:
- Outer lines placed at the k-th outermost point plus its expected gap k / (density × side length),
  with k chosen so that gap is 3 cm (so up to k−1 stray points beyond the edge are ignored).
  Median facet area error: gable −0.08%, hip +0.04%, cross-gable −0.07%. Median length error for
  every edge type ≤ 2 cm. Fallbacks 0; "other" sides: 2 (gable), 0, 0.
- If a polygon fails to close: drop the shortest crooked outer piece and retry (up to 4 times).
- Reject corners more than 2 m outside the traced outline (a spike). Good corners: up to 1.16 m
  on 200 synthetic facets, 1.44 m on the real houses; the delft_03 spike was 6.02 m.

### Real houses

Source: `scripts/stage3_topology.py` (RANSAC 0.12 m / 1 m², seed 0). "Fallback" facets keep the
traced outline (area from the outline, sides typed "other"). Lengths are totals; shared internal
edges counted once (average of the two facets' sides).

| house | facets | fallbacks | roof area m² polygons / hulls | ridge | hip | break | eave | rake | perimeter | other |
|---|---|---|---|---|---|---|---|---|---|---|
| delft_02 | 2 | 0 | 93.5 / 94.5 | | | | | | 55.5 | |
| delft_03 | 5 | 2 | 164.6 / 171.0 | 12.4 | | 2.3 | 24.9 | 17.7 | 24.3 | 19.5 |
| delft_05 | 4 | 2 | 92.8 / 102.9 | 7.4 | | | 14.7 | 23.1 | | 21.0 |
| delft_01 | 8 | 2 | 251.9 / 257.4 | | 12.3 | 23.5 | 39.7 | 2.4 | 40.7 | 35.2 |

- delft_03's gable: north half ridge 12.67, eave 12.67, rakes 4.45 / 4.40 m, 56.1 m² at 42.7°;
  south half ridge 12.22, eave 12.22, rakes 4.39 / 4.43 m, 53.9 m² at 44.1°. The four rakes agree
  within 6 cm. The shared ridge comes out 0.45 m different from each side: each facet ends it at
  its own rake corners. Needs reconciling into one edge.
- delft_03's flat extension reaches ~0.5 m under the gable's eave in plan. Probably real (the flat
  roof continues under the overhang and is partly seen at an angle), so plan areas overlap there.
- delft_01's mansard sides: eaves 11.3–12.0 m, hips ~4.0 m, breaks where they meet the flat top.
- Fallbacks are the small facets (2.9–5.0 m²: delft_03's end structure, delft_05's corner pieces)
  and delft_01's hump area and its 3.2 m² piece.
- delft_05's NW slope has two ridge sides (4.65 and 3.69 m), with the SE slope and a SE corner
  piece facing 18° apart: the twisting SE side from stage 2.

## Stage 4: measurement report

Source: `scripts/stage4_report.py` → `data/interim/<house>/report.json` (example:
`docs/examples/delft_03_report.json`). Each physical edge listed once; a shared edge's two versions
are averaged if they disagree by ≤ 0.5 m, otherwise their overlap is kept. 1 square = 100 ft².

Synthetic check (`tests/test_measure.py`, 3 seeds each): edge counts exactly gable 7, hip 9,
cross-gable 15; totals by type within 0.3 m; roof area within 2%.

| house | facets (fallback) | roof area m² | ft² | squares | ridge m | hip m | break m | eave m | rake m | perimeter m | unclassified m | max shared-edge disagreement m |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| delft_02 | 2 (0) | 93.5 | 1007 | 10.07 | | | | | | 55.46 | | |
| delft_03 | 5 (2) | 164.6 | 1772 | 17.72 | 12.45 | | 4.66 | 24.89 | 17.67 | 24.32 | 19.53 | 0.44 |
| delft_05 | 4 (2) | 92.8 | 999 | 9.99 | 7.60 | | | 14.73 | 23.08 | | 20.95 | 3.51 |
| delft_01 | 8 (2) | 251.9 | 2711 | 27.11 | | 12.03 | 23.45 | 39.67 | 2.39 | 40.69 | 35.24 | 0.58 |

- delft_05's ridge read 9.72 m when the two versions were averaged; its SE facet's ridge also runs
  along a corner piece, so the versions cover different stretches (3.51 m apart). Overlap: 7.60 m.
- "unclassified" = sides of fallback facets plus outer sides neither level nor up-slope. It is
  12–14% of the edge length on the three pitched houses: the small facets are the weak spot.
- Not yet validated against anything independent; that is stage 5.

## Stage 5: validation

### Footprint-based estimates

Source: `scripts/stage5_validate.py`. Footprint = BAG outline (walls). All estimates use the report's
pitches and plan-area shares, but not its areas, so they test consistency, not correctness.
A = footprint / cos(mean pitch); B = footprint × mean of 1/cos(pitch); C = B with the footprint widened
by the 98th-percentile overhang; D = B with footprint + overhang plan area ((points outside the
footprint) / density). Means weighted by plan area.

| house | report m² | A | B | C | D | overhang p98 / median m | overhang area m² |
|---|---|---|---|---|---|---|---|
| delft_01 | 227.8 | 183.0 (−19.7%) | 219.5 (−3.6%) | 240.5 (+5.6%) | 234.8 (+3.1%) | 0.23 / 0.09 | 11.9 |
| delft_02 | 93.5 | 92.3 (−1.3%) | 92.3 (−1.3%) | 102.1 (+9.2%) | 94.3 (+0.9%) | 0.22 / 0.06 | 2.0 |
| delft_03 | 164.6 | 136.7 (−17.0%) | 149.9 (−8.9%) | 183.3 (+11.4%) | 164.9 (+0.2%) | 0.57 / 0.17 | 12.0 |
| delft_05 | 92.5 | 91.0 (−1.6%) | 91.0 (−1.6%) | 104.8 (+13.3%) | 95.7 (+3.5%) | 0.33 / 0.07 | 4.4 |

- A vs B: 1/cos of the mean pitch understates a mixed-pitch roof (1/cos is convex). On delft_03
  (43° gable + flat extension) sec(mean pitch) = 1.141 vs mean sec = 1.251. On single-pitch roofs A = B.
- C overshoots by 6–13%: a uniform overhang taken from a high percentile assumes eaves all round;
  gable ends and flat roofs barely overhang (median outside-wall distance 6–17 cm).
- D, which counts the overhang's actual plan area, is within +0.2% to +3.5% of every report.
- This check found a report bug: on delft_01 two traced pieces of the curved top (20.8 + 3.2 m²) lay
  inside the flat top's polygon at the same height and were counted twice (report 251.9 m²; D said
  228.7). Fallback facets are now trimmed where they overlap a closed facet within 0.5 m in height;
  report 227.8 m². Remaining plan overlap (2.0–3.0 m² on delft_01/02/03) is a higher roof's edge over
  a lower roof: real, and kept.

### Synthetic checks: dormer

Roof: `tests/synthetic_roofs.dormer` (main gable 12 × 9 m at 40°, gable dormer 3 m wide on the south
slope, same pitch; dormer facets 3.82 m² each). 20 seeds, full pipeline (`build_roof` + report).
Facet area error = report area / true area − 1, median [min, max].

| version | planes found (of 4) | fallback facets | roof area | main south | dormer W | dormer E | valleys m (truth 4.93) | double-counted plan m² (median) |
|---|---|---|---|---|---|---|---|---|
| before | 4 in 19 runs, 3 in 1 | 8 | +5.7% | +12.0% | +17.3% [−40.5, +31.2] | +20.7% [−30.1, +32.5] | 6.46 | 6.3 |
| + facet holes | same | 8 | +0.8% | +0.6% | +17.3% | +20.7% | 6.46 | 1.5 |
| + edge sides must lie near the edge line | same | 8 | +0.6% | +0.6% | +4.5% [−40.5, +31.2] | +21.9% [−33.2, +32.5] | 6.28 | 1.3 |

Kept: both changes (all tests pass; gable/hip/cross-gable results unchanged).
- Holes: interior gaps of a facet's traced outline larger than 1 m² are subtracted. Without them the
  main slope's polygon covered the dormer.
- Beside a dormer's eave the main slope is just below, across the dormer's side wall, so the eave's
  outline "touched" the main slope and was taken for the valley; the valley line then ran past the
  eave corner. Outline points now count as on an internal edge only if they also lie within 0.3 m
  of its line.
- The run where RANSAC found 3 planes missed a dormer facet; its points joined the main slope, so
  there was no hole (main south +11.8%). A segmentation failure.

Tried and backed out (each improved the dormer and broke something else):
- Keeping short level/up-slope outer sides (≥ 0.3 m) instead of dropping all sides under 0.8 m:
  dormer W +0.3%, E still +20.7%, and a 0.78 m corner sliver on a hip became a "rake" (2 tests
  failed). Lengths of real short sides (dormer, 0.7–2.1 m) and of artefacts (plain roofs, up to
  2.35 m) overlap, so no length threshold separates them.
- Alpha radius sweep 0.35–0.75 m: dormer facet |error| median 22.5–25.9% at every value. Not the cause.
- A second candidate polygon per facet, the intersection of half-planes (its internal lines plus
  level/up-slope lines at the points' extremes), chosen against the run-based polygon by IoU with the
  traced outline: dormer facets −5.2% / −2.1% median, but fallbacks 8 → 13, valleys 4.93 → 2.62 m,
  a −6% hip facet, 4 tests failed. Likely flaw: the traced outline has rounded corners, so IoU against
  it favours polygons with clipped corners.

Open: dormer facets (≈ 4 m², sides 1.2–2 m) are not measured reliably: dormer E +21.9% median,
spread −33% to +32%, 8 of 40 fall back. This matters for Bilthoven, where dormers are common.

Real houses after these changes: delft_01 224.2 m² (was 227.8; the flat top now has holes where the
curved parts are), delft_02 93.5, delft_03 164.4, delft_05 92.9. Against footprint estimate D: +4.3%,
+0.9%, +0.4%, +3.0%.

### Area uncertainty by facet size

Source: full pipeline + report on the synthetic gable, hip, cross-gable and dormer roofs, 20 seeds
each. |area error| by true facet area:

| true facet area | facets | median | 90th percentile | max | fallbacks |
|---|---|---|---|---|---|
| < 10 m² | 39 | 23.0% | 31.5% | 40.5% | 8 |
| 10–25 m² | 80 | 0.5% | 1.2% | 1.6% | 0 |
| ≥ 25 m² | 160 | 0.2% | 0.7% | 11.8% | 0 |

No synthetic facet lies between 3.8 and 17.1 m², so the 10 m² boundary is a choice within that
range. The ≥ 25 m² maximum (11.8%) is the run where RANSAC missed a dormer facet. These are ideal
roofs, so real errors are at least this large.

The report now gives each facet a confidence (low if < 10 m² or traced) and its 90th-percentile
error, and a conservative total bound (facet errors added as if all in the same direction):

| house | roof area m² | ± p90 bound m² | in low-confidence facets m² |
|---|---|---|---|
| delft_01 | 224.2 | 6.0 | 14.2 (6%) |
| delft_02 | 93.5 | 0.7 | 0 |
| delft_03 | 164.4 | 4.1 | 9.5 (6%) |
| delft_05 | 92.9 | 3.0 | 7.6 (8%) |
