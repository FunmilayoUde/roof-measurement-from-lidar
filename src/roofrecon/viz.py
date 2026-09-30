"""Rendering: a static matplotlib snapshot (always works, diffable) and an Open3D viewer."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def save_snapshot(xyz: np.ndarray, out_path: str | Path, footprint_xy: np.ndarray | None = None, title: str = "") -> None:
    """Plan view coloured by height plus two elevations. Enough to see ridges, pitch and outliers."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    x, y, z = xyz.T
    sc = axes[0].scatter(x, y, c=z, s=2, cmap="viridis")
    if footprint_xy is not None:
        axes[0].plot(footprint_xy[:, 0], footprint_xy[:, 1], "r-", lw=1, label="footprint")
        axes[0].legend(loc="upper right")
    axes[0].set_aspect("equal")
    axes[0].set_title("plan (colour = z)")
    fig.colorbar(sc, ax=axes[0], label="z [m]")
    for ax, (h, label) in zip(axes[1:], [(x, "x"), (y, "y")]):
        ax.scatter(h, z, c=z, s=1, cmap="viridis")
        ax.set_aspect("equal")
        ax.set_xlabel(f"{label} [m]")
        ax.set_ylabel("z [m]")
        ax.set_title(f"elevation along {label}")
    fig.suptitle(title)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def to_open3d(xyz: np.ndarray):
    """Open3D point cloud coloured by height."""
    import matplotlib
    import open3d as o3d

    pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(xyz))
    z = xyz[:, 2]
    t = (z - z.min()) / max(np.ptp(z), 1e-9)
    pcd.colors = o3d.utility.Vector3dVector(matplotlib.colormaps["viridis"](t)[:, :3])
    return pcd


def write_ply(xyz: np.ndarray, path: str | Path) -> None:
    import open3d as o3d

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not o3d.io.write_point_cloud(str(path), to_open3d(xyz), write_ascii=False):
        raise OSError(f"Open3D failed to write {path}")
