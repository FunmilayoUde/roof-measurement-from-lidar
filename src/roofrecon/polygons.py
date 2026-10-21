"""Facet polygons: each facet's outline as straight sides, with corners where the sides' lines meet.

Each facet is handled in its own 2D plane coordinates (u horizontal, v up the slope; see
planes.plane_axes). A side is either
  - internal: it runs along an edge shared with another facet, and its line is that exact
    plane-plane intersection (topology.internal_edges), or
  - outer: on the roof's outline. Its line is fitted to the facet's boundary points and
    snapped to horizontal (u) or straight up-slope (v) when close, because eaves are level
    and gable rakes run straight up the slope.
Corners are intersections of consecutive side lines, so they don't depend on where the
outermost points happen to fall.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import shapely
from scipy.spatial import cKDTree
from shapely.geometry import Polygon

from roofrecon.planes import Plane, alpha_shape, plane_axes, to_plane_coords
from roofrecon.segment import Segmentation
from roofrecon.topology import Edge, plane_intersection

FLAT_PITCH_SIN = np.sin(np.radians(5.0))


@dataclass
class Side:
    kind: str  # internal edge kind (ridge|hip|valley|break) or outer kind (eave|rake|top|perimeter|other)
    neighbour: int | None  # the other facet for internal sides
    p0: np.ndarray  # 3D
    p1: np.ndarray

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))


@dataclass
class FacetPolygon:
    facet: int
    plane: Plane
    vertices: np.ndarray  # (K, 3), in order, not repeated at the end
    sides: list[Side] = field(default_factory=list)
    fallback: bool = False  # True if line intersections failed and the traced outline was used

    @property
    def area(self) -> float:
        """True (sloped) area: shoelace area in plane coordinates."""
        uv = to_plane_coords(self.vertices, self.plane)
        return float(Polygon(uv).area)

    @property
    def plan_area(self) -> float:
        return float(Polygon(self.vertices[:, :2]).area)


def _line_from_3d(point, direction, plane: Plane):
    """A 3D line lying in the plane, as a 2D (point, unit direction) in plane coordinates."""
    u, v = plane_axes(plane.normal)
    p = to_plane_coords(point[None, :], plane)[0]
    d = np.array([direction @ u, direction @ v])
    return p, d / np.linalg.norm(d)


def _intersect_2d(l1, l2, min_sin: float = np.sin(np.radians(10))):
    (p, d), (q, e) = l1, l2
    cross = d[0] * e[1] - d[1] * e[0]
    if abs(cross) < min_sin:
        return None
    t = ((q - p)[0] * e[1] - (q - p)[1] * e[0]) / cross
    return p + t * d


def _to_3d(uv: np.ndarray, plane: Plane) -> np.ndarray:
    u, v = plane_axes(plane.normal)
    return plane.centroid + uv[:, :1] * u + uv[:, 1:2] * v


def _fit_outer_line(pts: np.ndarray, snap_deg: float):
    """Total-least-squares line through 2D points, snapped to u (horizontal) or v (up-slope)
    if within snap_deg."""
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c)
    d = vt[0]
    for axis in (np.array([1.0, 0.0]), np.array([0.0, 1.0])):
        if abs(d @ axis) > np.cos(np.radians(snap_deg)):
            d = axis
            break
    return c, d


def _push_to_edge(line, uv: np.ndarray, centre: np.ndarray, run_pts: np.ndarray, rank: int = 2):
    """Move an outer line outward to the facet's outermost points along it.

    The line was fitted through the outline's vertices, which are boundary points: they sit
    inside the true edge by a few cm, since random samples rarely land right on it (on
    synthetic gables, area read 2.7% low). Among the facet's points alongside this side,
    the `rank`-th farthest out is a near-unbiased estimate of the edge; using the 2nd
    rather than the 1st keeps a single stray point from pushing the edge out.
    """
    p, d = line
    normal = np.array([-d[1], d[0]])
    if (centre - p) @ normal > 0:  # make the normal point outward, away from the facet
        normal = -normal
    t = (run_pts - p) @ d
    alongside = ((uv - p) @ d >= t.min()) & ((uv - p) @ d <= t.max())
    out = np.sort((uv[alongside] - p) @ normal)[::-1]
    if len(out) < rank:
        return line
    return p + out[rank - 1] * normal, d


def _outer_kind(line, facet_uv_mean, flat: bool) -> str:
    """eave: level side below the facet; top: level side above it; rake: side running up the
    slope; perimeter: any outer side of a flat facet; other: anything else."""
    p, d = line
    if flat:
        return "perimeter"
    if abs(d[1]) < 1e-9:  # horizontal in 3D
        return "eave" if p[1] < facet_uv_mean[1] else "top"
    if abs(d[0]) < 1e-9:
        return "rake"
    return "other"


def facet_polygon(
    xyz: np.ndarray,
    seg: Segmentation,
    facet: int,
    edges: list[Edge],
    alpha_radius: float = 0.75,
    simplify_tol: float = 0.25,
    touch_radius: float = 0.4,
    snap_deg: float = 10.0,
    min_side: float = 0.8,
    trim: float = 0.4,
) -> FacetPolygon:
    plane = seg.planes[facet]
    pts = xyz[seg.labels == facet]
    uv = to_plane_coords(pts, plane)
    flat = np.linalg.norm(plane.normal[:2]) < FLAT_PITCH_SIN

    outline = alpha_shape(uv, alpha_radius)
    if outline.geom_type == "MultiPolygon":
        outline = max(outline.geoms, key=lambda g: g.area)
    ring = np.asarray(outline.exterior.coords)[:-1]
    simple = np.asarray(shapely.simplify(Polygon(ring), simplify_tol).exterior.coords)[:-1]

    # Internal lines available to this facet, keyed by the neighbouring facet.
    internal = {}
    for e in edges:
        if facet not in e.facets:
            continue
        other = e.facets[1] if e.facets[0] == facet else e.facets[0]
        line = plane_intersection(seg.planes[e.facets[0]], seg.planes[e.facets[1]])
        if line is not None:
            internal[other] = (e, _line_from_3d(*line, plane))

    # Label each outline vertex by the neighbouring facet it touches. This is topological:
    # a stretch of outline that borders the wing's points is the edge shared with the wing,
    # however ragged the labels make it look. (Matching outline segments to edge lines by
    # distance and angle failed exactly there.)
    ring3d = _to_3d(ring, plane)
    near = np.full(len(ring), -1)
    best = np.full(len(ring), np.inf)
    for other in internal:
        dist, _ = cKDTree(xyz[seg.labels == other]).query(ring3d)
        closer = (dist < touch_radius) & (dist < best)
        near[closer], best[closer] = other, dist[closer]
    near = _majority_smooth(near, window=5)

    # Split the ring into runs of equal label; split outer stretches further into straight
    # pieces where the outline turns (an eave running into a rake, say).
    runs = []
    for label, pts in _label_runs(ring, near):
        if label >= 0:
            e, line = internal[label]
            runs.append({"label": ("internal", e, label, line), "pts": pts})
            continue
        if np.allclose(pts[0], pts[-1]):
            # The whole outline, touching no other facet: simplify it as a ring, then rotate
            # it to start at a kept corner so the pieces can be cut in order.
            loop = pts[:-1]
            corners = np.asarray(shapely.simplify(Polygon(loop), simplify_tol).exterior.coords)[:-1]
            idx = sorted({int(np.argmin(np.linalg.norm(loop - q, axis=1))) for q in corners})
            pts = np.vstack([loop[idx[0]:], loop[:idx[0] + 1]])
            idx = [i - idx[0] for i in idx] + [len(pts) - 1]
        else:
            line_simple = np.asarray(shapely.simplify(shapely.linestrings(pts), simplify_tol).coords)
            idx = [int(np.argmin(np.linalg.norm(pts - q, axis=1))) for q in line_simple]
        for i0, i1 in zip(idx[:-1], idx[1:]):
            if i1 > i0:
                runs.append({"label": ("outer",), "pts": pts[i0:i1 + 1]})

    # Merge consecutive runs with the same internal edge, or consecutive outer runs that
    # are nearly collinear.
    merged = []
    for r in runs:
        if merged and _same_label(merged[-1], r, snap_deg):
            merged[-1]["pts"] = np.vstack([merged[-1]["pts"], r["pts"][1:]])
            continue
        merged.append(r)
    if len(merged) > 1 and _same_label(merged[-1], merged[0], snap_deg):
        merged[0]["pts"] = np.vstack([merged[-1]["pts"], merged[0]["pts"][1:]])
        merged.pop()

    # The traced outline rounds every corner, and simplification turns each rounding into a
    # short diagonal side. Drop short runs so their neighbours meet in a sharp corner, then
    # merge again, since dropping one can leave two collinear outer runs side by side.
    merged = _drop_short_runs(merged, min_side, snap_deg)

    # A line per run. Consecutive runs whose lines come out parallel can't make a corner:
    # two outer runs are merged and refitted; an outer run beside a parallel internal one
    # is dropped (the internal line is exact).
    centre = uv.mean(axis=0)

    def line_of(r):
        if r["label"][0] == "internal":
            return r["label"][3]
        line = _fit_outer_line(_trim_ends(r["pts"], trim), snap_deg)
        return _push_to_edge(line, uv, centre, r["pts"])

    changed = True
    while changed and len(merged) > 3:
        changed = False
        for i in range(len(merged)):
            j = (i + 1) % len(merged)
            if _intersect_2d(line_of(merged[i]), line_of(merged[j])) is not None:
                continue
            a_int, b_int = merged[i]["label"][0] == "internal", merged[j]["label"][0] == "internal"
            if a_int and b_int:
                continue  # two parallel internal lines: leave it to the fallback check
            if not a_int and not b_int:
                merged[i]["pts"] = np.vstack([merged[i]["pts"], merged[j]["pts"]])
                merged.pop(j)
            else:
                merged.pop(j if a_int else i)
            changed = True
            break

    lines, kinds, neighbours = [], [], []
    for r in merged:
        line = line_of(r)
        lines.append(line)
        if r["label"][0] == "internal":
            kinds.append(r["label"][1].kind)
            neighbours.append(r["label"][2])
        else:
            kinds.append(_outer_kind(line, uv.mean(axis=0), flat))
            neighbours.append(None)

    # Corners: intersection of each run's line with the next one's.
    corners = []
    ok = len(lines) >= 3
    for i in range(len(lines)):
        q = _intersect_2d(lines[i], lines[(i + 1) % len(lines)]) if ok else None
        if q is None:
            ok = False
            break
        corners.append(q)
    if ok:
        poly2d = np.array(corners)
        candidate = Polygon(poly2d)
        # Reject corner sets that went wild (self-intersecting, or far from the traced outline).
        if not candidate.is_valid or abs(candidate.area - outline.area) > 0.5 * outline.area:
            ok = False
    if not ok:
        poly2d = simple
        verts = _to_3d(poly2d, plane)
        sides = [Side("other", None, verts[i], verts[(i + 1) % len(verts)]) for i in range(len(verts))]
        return FacetPolygon(facet, plane, verts, sides, fallback=True)

    verts = _to_3d(poly2d, plane)
    # Side i runs from corner i-1 to corner i (corner i ends run i).
    sides = [Side(kinds[i], neighbours[i], verts[i - 1], verts[i]) for i in range(len(lines))]
    return FacetPolygon(facet, plane, verts, sides)


def _majority_smooth(labels: np.ndarray, window: int) -> np.ndarray:
    """Circular majority filter: removes single-vertex flicker between labels."""
    n, h = len(labels), window // 2
    out = labels.copy()
    for i in range(n):
        vals, counts = np.unique(labels[[(i + k) % n for k in range(-h, h + 1)]], return_counts=True)
        out[i] = vals[np.argmax(counts)]
    return out


def _label_runs(ring: np.ndarray, labels: np.ndarray) -> list[tuple[int, np.ndarray]]:
    """Consecutive vertices with the same label, going round the ring once. Each run includes
    the first vertex of the next run, so consecutive runs share their joining vertex."""
    n = len(ring)
    changes = np.flatnonzero(labels != np.roll(labels, 1))
    if len(changes) == 0:
        return [(int(labels[0]), np.vstack([ring, ring[:1]]))]
    out = []
    for a, b in zip(changes, np.append(changes[1:], changes[0] + n)):
        idx = np.arange(a, b + 1) % n
        out.append((int(labels[a]), ring[idx]))
    return out


def _run_length(r) -> float:
    return float(np.linalg.norm(r["pts"][-1] - r["pts"][0]))


def _same_label(a, b, deg) -> bool:
    if a["label"][0] == "internal" and b["label"][0] == "internal":
        return a["label"][1] is b["label"][1]
    if a["label"][0] == "outer" and b["label"][0] == "outer":
        return _collinear(a["pts"], b["pts"], deg)
    return False


def _is_artefact(runs: list[dict], i: int, min_side: float, deg: float, max_artefact: float) -> bool:
    """Is run i a rounding of the traced outline rather than a real side?
    - any run shorter than min_side;
    - an outer run shorter than max_artefact that is neither level nor straight up-slope
      (a corner cut off diagonally);
    - an outer run shorter than max_artefact between two internal runs (the tip of a notch,
      where the two internal lines should meet at a three-plane corner)."""
    r = runs[i]
    length = _run_length(r)
    if length < min_side:
        return True
    if r["label"][0] != "outer" or length >= max_artefact:
        return False
    d = _fit_outer_line(r["pts"], 0.0)[1]
    snapped = max(abs(d[0]), abs(d[1])) > np.cos(np.radians(deg))
    prev_int = runs[i - 1]["label"][0] == "internal"
    next_int = runs[(i + 1) % len(runs)]["label"][0] == "internal"
    return not snapped or (prev_int and next_int)


def _drop_short_runs(runs: list[dict], min_side: float, deg: float, max_artefact: float = 2.0) -> list[dict]:
    runs = list(runs)
    while len(runs) > 3:
        candidates = [i for i in range(len(runs)) if _is_artefact(runs, i, min_side, deg, max_artefact)]
        if not candidates:
            break
        i = min(candidates, key=lambda i: _run_length(runs[i]))
        runs.pop(i)
        # Re-merge the neighbours that are now adjacent, if they belong together.
        j = i % len(runs)
        k = (j - 1) % len(runs)
        if len(runs) > 3 and _same_label(runs[k], runs[j], deg):
            runs[k]["pts"] = np.vstack([runs[k]["pts"], runs[j]["pts"]])
            runs.pop(j)
    return runs


def _trim_ends(pts: np.ndarray, trim: float) -> np.ndarray:
    """Drop points within `trim` of either end of a run (where the outline rounds a corner),
    keeping at least the middle few points."""
    if len(pts) < 5:
        return pts
    along = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    keep = (along > trim) & (along < along[-1] - trim)
    return pts[keep] if keep.sum() >= 3 else pts[len(pts) // 4: -(len(pts) // 4) or None]


def _collinear(pts_a: np.ndarray, pts_b: np.ndarray, deg: float) -> bool:
    da = _fit_outer_line(pts_a, 0.0)[1] if len(pts_a) >= 2 else None
    db = _fit_outer_line(pts_b, 0.0)[1] if len(pts_b) >= 2 else None
    if da is None or db is None:
        return True
    return abs(da @ db) > np.cos(np.radians(deg))


def roof_polygons(xyz: np.ndarray, seg: Segmentation, edges: list[Edge], **kw) -> list[FacetPolygon]:
    return [facet_polygon(xyz, seg, k, edges, **kw) for k in range(len(seg.planes))]
