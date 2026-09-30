"""Building footprints from OpenStreetMap, reprojected into the point cloud's CRS."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

from shapely.geometry import Polygon, shape
from shapely.ops import transform

OVERPASS_URL = "https://overpass-api.de/api/interpreter"


def fetch_osm_building(lat: float, lon: float, out_path: str | Path) -> dict:
    """Fetch the OSM building way that contains (lat, lon) and save it as GeoJSON (EPSG:4326).

    `is_in` + `pivot` returns areas enclosing the point, which is what we want: a point
    inside a house is usually far from its outline, so a distance-based `around` query
    would miss it.
    """
    query = f"""
    [out:json][timeout:25];
    is_in({lat},{lon})->.a;
    way(pivot.a)["building"];
    out geom tags;
    """
    req = urllib.request.Request(
        OVERPASS_URL,
        data=urllib.parse.urlencode({"data": query}).encode(),
        headers={"User-Agent": "roofrecon/0.1 (learning project)"},
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.load(resp)

    ways = [e for e in payload.get("elements", []) if e.get("type") == "way"]
    if not ways:
        raise ValueError(f"No OSM building way contains ({lat}, {lon})")
    if len(ways) > 1:
        # Nested building parts; the smallest enclosing way is the most specific.
        ways.sort(key=lambda w: Polygon([(p["lon"], p["lat"]) for p in w["geometry"]]).area)
    way = ways[0]

    feature = {
        "type": "Feature",
        "properties": {"osm_id": way["id"], **way.get("tags", {})},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[p["lon"], p["lat"]] for p in way["geometry"]]],
        },
    }
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(feature, indent=2))
    return feature


def load_footprint(path: str | Path, src_crs: str = "EPSG:4326", dst_crs: str | None = None) -> Polygon:
    """Load a single-polygon GeoJSON Feature and reproject it (x, y order) into dst_crs."""
    feature = json.loads(Path(path).read_text())
    geom = shape(feature["geometry"] if "geometry" in feature else feature)
    if dst_crs is None or dst_crs == src_crs:
        return geom
    from pyproj import Transformer

    # always_xy: EPSG:4326 is formally (lat, lon); forcing (lon, lat) avoids axis swaps.
    t = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    return transform(t.transform, geom)
