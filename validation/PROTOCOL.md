# Hand-measurement protocol

The check that the pipeline can't grade itself on: you measure roof features by hand, then
`scripts/stage5_validate.py` compares them with each house's `report.json` and reports the error
distribution.

## Rule 1: measure before you look

Don't open `report.json`, `report_plan.png`, or the stage 2–4 previews for a house until you've
finished measuring it. Knowing the answer pulls a measurement towards it. Matching your
measurements to report edges (step 4) comes after.

## What to measure, per house

Aim for 3–4 houses done carefully rather than many done quickly. Per house:

- **Every eave and rake you can see clearly** (lengths).
- **The ridge(s)**, and hips or valleys if there are any (lengths).
- **The pitch of each main slope.**
- Optionally, **plan lengths from the aerial photo** (method 2) for the same edges.

## Method 1: CloudCompare, on the point cloud

Independent of the algorithm, though not of the data (same LiDAR points).

1. Open `data/interim/<house>/roof.ply`. Its coordinates are local (recentered in stage 1); if
   CloudCompare offers a global shift, decline it, so the coordinates you read match the report's.
2. Use the point-picking tool in its two-point (distance) mode.
3. **Edge length:** pick a point at each end of the edge, on the edge itself. Roof points don't sit
   exactly on edges; pick the outermost points that are clearly on the line. Record both points'
   x, y, z. The script computes the length.
4. **Pitch:** pick two points on the same slope, one near the top and one near the bottom, straight
   down the slope (not diagonally). Record both points. The script computes rise over run.
5. Pick each pair twice, independently, and record both rows. The difference between your own
   repeats tells us how precise a hand measurement is, which bounds how small a pipeline error we
   can detect.

## Method 2: aerial photo (plan lengths only)

Independent of the LiDAR entirely: a different sensor and a different date.

- Open the house in the PDOK viewer (`app.pdok.nl`, layer "Luchtfoto actueel"), or in QGIS with
  the PDOK aerial photo WMTS layer, and measure along roof edges.
- Measure **lengths, not positions**. The photo isn't a true orthophoto, so tall parts lean away
  from the image centre and roof outlines are shifted by up to ~0.5 m, but an edge's length is
  mostly preserved.
- A sloped edge (rake, hip) shows its plan length here, shorter than its true length. The script
  compares plan length with plan length.

## Recording: `validation/hand_measurements.csv`

One row per measurement. Columns:

| column | meaning |
|---|---|
| `house_id` | e.g. `delft_03` |
| `kind` | `edge` (two picked points), `slope` (two picked points, top then bottom), or `plan_length` (photo) |
| `feature` | your own description, e.g. "north eave of main gable" |
| `x0,y0,z0,x1,y1,z1` | the two picked points (local coordinates); empty for `plan_length` |
| `value` | for `plan_length`: the length in m; empty otherwise |
| `report_id` | filled in step 4: the report's edge id (edges) or facet id (slopes). For `edge` and `slope` rows it can stay empty and the script matches by position. |
| `method` | `cloudcompare` or `photo` |
| `repeat` | 1 or 2 |
| `notes` | anything odd: trees, a hard-to-see corner, uncertain end |

Example rows:

```
house_id,kind,feature,x0,y0,z0,x1,y1,z1,value,report_id,method,repeat,notes
delft_03,edge,north eave of main gable,-6.9,2.7,5.6,4.9,7.1,5.6,,,cloudcompare,1,
delft_03,slope,north slope of main gable,0.1,3.0,8.4,-0.9,5.6,5.9,,,cloudcompare,1,
delft_03,plan_length,north eave of main gable,,,,,,,12.6,7,photo,1,measured on 2024 photo
```

## Step 4: matching

After measuring a house, run `uv run python scripts/stage4_report.py --house-id <house>` and open
`data/interim/<house>/report_plan.png`, which numbers every edge and facet. Fill in `report_id` for
`plan_length` rows (required) and for any `edge` or `slope` rows the script matched wrongly.

## Then

```sh
uv run python scripts/stage5_validate.py
```

prints each measurement next to the report's value, and the error distribution (median, spread,
worst case) per kind. The numbers go in RESULTS.md.
