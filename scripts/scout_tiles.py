"""Score GeoTiles AHN sub-tiles around a place by how many useful detached houses they hold.

    uv run python scripts/scout_tiles.py --place wassenaar 52.143 4.401 --place bilthoven 52.130 5.205

AHN sheets are 5 km x 6.25 km in EPSG:28992, aligned to multiples of 5000 (x) and 6250 (y);
GeoTiles splits each into 25 sub-tiles of 1 km x 1.25 km, numbered row by row from the
north-west (01) to the south-east (25). This finds the 3x3 sub-tiles around each place and
counts candidate detached houses (as find_candidates.py), and how many have L/T-shaped plans
(concavity > 0.1), which is where valleys come from. Sheet names (e.g. 37EN2) aren't derived
here; look them up on GeoTiles from the printed coordinates.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer
from shapely.geometry import box, shape
from shapely.strtree import STRtree

sys.path.insert(0, str(Path(__file__).parent))
from find_candidates import bag_buildings_in_bbox  # noqa: E402

SHEET_W, SHEET_H, SUB_W, SUB_H = 5000, 6250, 1000, 1250


def subtile(x: float, y: float) -> tuple[float, float, int]:
    """Sheet origin (south-west corner) and sub-tile number (1-25) containing a point."""
    x0, y0 = SHEET_W * np.floor(x / SHEET_W), SHEET_H * np.floor(y / SHEET_H)
    col = int((x - x0) // SUB_W) + 1
    row = int((y0 + SHEET_H - y) // SUB_H) + 1
    return x0, y0, (row - 1) * 5 + col


def score(bounds) -> dict:
    tile = box(*bounds)
    feats = bag_buildings_in_bbox(bounds)
    polys = [shape(f["geometry"]) for f in feats]
    tree = STRtree(polys)
    n = concave = 0
    for i, (f, p) in enumerate(zip(feats, polys)):
        pr = f["properties"]
        if not (
            "woonfunctie" in (pr.get("gebruiksdoel") or "")
            and pr.get("aantal_verblijfsobjecten") == 1
            and 90 <= p.area <= 400
            and tile.contains(p.buffer(5))
        ):
            continue
        if any(polys[j].distance(p) < 3.0 for j in tree.query(p.buffer(3.0)) if j != i):
            continue
        n += 1
        concave += (1 - p.area / p.convex_hull.area) > 0.1
    return {"buildings": len(feats), "candidates": n, "l_or_t_shaped": concave}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--place", nargs=3, action="append", metavar=("NAME", "LAT", "LON"), required=True)
    args = ap.parse_args()
    to_rd = Transformer.from_crs("EPSG:4326", "EPSG:28992", always_xy=True)
    for name, lat, lon in args.place:
        cx, cy = to_rd.transform(float(lon), float(lat))
        print(f"{name} (RD {cx:.0f}, {cy:.0f})")
        for dy in (SUB_H, 0, -SUB_H):
            for dx in (-SUB_W, 0, SUB_W):
                x, y = cx + dx, cy + dy
                x0, y0, k = subtile(x, y)
                col, row = (k - 1) % 5, (k - 1) // 5
                bx0, by1 = x0 + col * SUB_W, y0 + SHEET_H - row * SUB_H
                bounds = (bx0, by1 - SUB_H, bx0 + SUB_W, by1)
                s = score(bounds)
                print(
                    f"  sheet SW ({x0:.0f}, {y0:.0f}) sub-tile {k:02d}  bbox {tuple(int(b) for b in bounds)}: "
                    f"{s['buildings']:>4} buildings, {s['candidates']:>3} detached candidates, "
                    f"{s['l_or_t_shaped']:>3} L/T-shaped"
                )


if __name__ == "__main__":
    main()
