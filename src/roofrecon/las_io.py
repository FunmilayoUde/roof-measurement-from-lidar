"""Reading LAS/LAZ tiles. laspy is imported lazily so the geometry core stays numpy-only."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class PointData:
    xyz: np.ndarray  # (N, 3) float64, in the tile's CRS
    classification: np.ndarray  # (N,) uint8, ASPRS codes
    return_number: np.ndarray  # (N,) uint8
    number_of_returns: np.ndarray  # (N,) uint8
    header: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.xyz)

    def subset(self, mask: np.ndarray) -> PointData:
        return PointData(
            self.xyz[mask],
            self.classification[mask],
            self.return_number[mask],
            self.number_of_returns[mask],
            self.header,
        )


def header_summary(path: str | Path) -> dict:
    """Header facts that matter for geometry: extent, count, and the integer quantization step."""
    import laspy

    with laspy.open(path) as f:
        h = f.header
        crs = None
        try:
            parsed = h.parse_crs()
            crs = parsed.to_string() if parsed is not None else None
        except Exception:  # CRS parsing needs pyproj and a well-formed VLR; neither is guaranteed
            pass
        return {
            "path": str(path),
            "point_format": h.point_format.id,
            "version": str(h.version),
            "point_count": int(h.point_count),
            "scales": [float(s) for s in h.scales],
            "offsets": [float(o) for o in h.offsets],
            "mins": [float(v) for v in h.mins],
            "maxs": [float(v) for v in h.maxs],
            "crs": crs,
        }


def read_points(
    path: str | Path,
    bbox: tuple[float, float, float, float] | None = None,
    chunk_size: int = 2_000_000,
) -> PointData:
    """Stream a tile in chunks, keeping only points inside bbox = (xmin, ymin, xmax, ymax).

    A 1 km tile is tens of millions of points; one house is a few thousand. Chunked reading
    keeps memory flat instead of materialising the whole tile.
    """
    import laspy

    xyz_parts, cls_parts, rn_parts, nr_parts = [], [], [], []
    with laspy.open(path) as f:
        for chunk in f.chunk_iterator(chunk_size):
            x = np.asarray(chunk.x, dtype=np.float64)
            y = np.asarray(chunk.y, dtype=np.float64)
            if bbox is not None:
                xmin, ymin, xmax, ymax = bbox
                m = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
            else:
                m = np.ones(len(x), dtype=bool)
            if not m.any():
                continue
            z = np.asarray(chunk.z, dtype=np.float64)
            xyz_parts.append(np.column_stack([x[m], y[m], z[m]]))
            cls_parts.append(np.asarray(chunk.classification, dtype=np.uint8)[m])
            rn_parts.append(np.asarray(chunk.return_number, dtype=np.uint8)[m])
            nr_parts.append(np.asarray(chunk.number_of_returns, dtype=np.uint8)[m])

    if not xyz_parts:
        empty = np.empty(0, dtype=np.uint8)
        return PointData(np.empty((0, 3)), empty, empty, empty, header_summary(path))
    return PointData(
        np.concatenate(xyz_parts),
        np.concatenate(cls_parts),
        np.concatenate(rn_parts),
        np.concatenate(nr_parts),
        header_summary(path),
    )
