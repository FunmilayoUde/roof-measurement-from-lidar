"""List detached-house candidates inside a LAZ tile, from Dutch BAG footprints.

    uv run python scripts/find_candidates.py --laz data/raw/<tile>.LAZ

A candidate is a residential BAG building with one dwelling, a house-sized footprint,
no other building within --min-gap metres, and completely inside the tile. The gap is
measured between wall outlines, so it has to exceed two eave overhangs for the roofs
not to touch. Writes data/interim/candidates_<tile>.csv, nearest-first by gap.
"""

from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

from shapely.geometry import box, shape
from shapely.strtree import STRtree

from roofrecon.footprint import BAG_WFS_URL
from roofrecon.las_io import header_summary

LOCATIESERVER = "https://api.pdok.nl/bzk/locatieserver/search/v3_1/reverse"


def bag_buildings_in_bbox(bbox: tuple[float, float, float, float], page: int = 1000) -> list[dict]:
    """All BAG 'pand' features in an EPSG:28992 bbox. The WFS caps each response, so page."""
    feats, start = [], 0
    while True:
        params = {
            "service": "WFS", "version": "2.0.0", "request": "GetFeature", "typeNames": "bag:pand",
            "outputFormat": "application/json", "srsName": "EPSG:28992",
            "bbox": ",".join(map(str, bbox)) + ",EPSG:28992",
            "count": page, "startIndex": start, "sortBy": "identificatie",
        }
        with urllib.request.urlopen(f"{BAG_WFS_URL}?{urllib.parse.urlencode(params)}", timeout=120) as r:
            batch = json.load(r)["features"]
        feats += batch
        if len(batch) < page:
            return feats
        start += page


def nearest_address(x: float, y: float) -> str:
    from pyproj import Transformer

    lon, lat = Transformer.from_crs("EPSG:28992", "EPSG:4326", always_xy=True).transform(x, y)
    q = urllib.parse.urlencode({"lat": lat, "lon": lon, "type": "adres", "rows": 1})
    with urllib.request.urlopen(f"{LOCATIESERVER}?{q}", timeout=30) as r:
        docs = json.load(r)["response"]["docs"]
    return docs[0]["weergavenaam"] if docs else ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--laz", required=True, type=Path)
    ap.add_argument("--min-area", type=float, default=70.0)
    ap.add_argument("--max-area", type=float, default=300.0)
    ap.add_argument("--min-gap", type=float, default=3.0)
    args = ap.parse_args()

    h = header_summary(args.laz)
    tile = box(h["mins"][0], h["mins"][1], h["maxs"][0], h["maxs"][1])
    feats = bag_buildings_in_bbox(tile.bounds)
    polys = [shape(f["geometry"]) for f in feats]
    tree = STRtree(polys)
    print(f"{len(feats)} BAG buildings in tile")

    rows = []
    for i, (f, p) in enumerate(zip(feats, polys)):
        pr = f["properties"]
        if not (
            "woonfunctie" in (pr.get("gebruiksdoel") or "")
            and pr.get("aantal_verblijfsobjecten") == 1
            and args.min_area <= p.area <= args.max_area
            and tile.contains(p.buffer(5))
        ):
            continue
        gaps = sorted(polys[j].distance(p) for j in tree.query(p.buffer(30)) if j != i)
        gap = gaps[0] if gaps else 30.0
        if gap < args.min_gap:
            continue
        rows.append({
            "bag_id": pr["identificatie"], "area_m2": round(p.area, 1), "vertices": len(p.exterior.coords) - 1,
            "gap_m": round(gap, 1), "bouwjaar": pr.get("bouwjaar"),
            "x": round(p.centroid.x, 1), "y": round(p.centroid.y, 1),
        })
    rows.sort(key=lambda r: -r["gap_m"])
    for r in rows:
        r["address"] = nearest_address(r["x"], r["y"])

    out = Path("data/interim") / f"candidates_{args.laz.stem}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["bag_id"])
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} candidates -> {out}")
    for r in rows:
        print(f"  {r['bag_id']}  {r['area_m2']:>6} m²  {r['vertices']:>3} vtx  gap {r['gap_m']:>5} m  {r['bouwjaar']}  {r['address']}")


if __name__ == "__main__":
    main()
