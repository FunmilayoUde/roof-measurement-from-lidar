"""Rendering: a static matplotlib snapshot (always works, diffable) and an Open3D viewer."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def save_snapshot(
    xyz: np.ndarray, out_path: str | Path, footprint_xy: np.ndarray | None = None, title: str = ""
) -> None:
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


def segment_colors(labels: np.ndarray) -> np.ndarray:
    """RGB per point: a distinct colour per plane, grey for unassigned (-1)."""
    import matplotlib

    palette = matplotlib.colormaps["tab20"].colors
    colors = np.array([palette[k % len(palette)] if k >= 0 else (0.75, 0.75, 0.75) for k in labels])
    return colors


def save_segment_snapshot(xyz: np.ndarray, labels: np.ndarray, out_path: str | Path, title: str = "") -> None:
    """Plan and two elevations coloured by plane, with plane numbers at their centres in plan."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = segment_colors(labels)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    x, y, z = xyz.T
    for ax, (h, v) in zip(axes, [(x, y), (x, z), (y, z)]):
        ax.scatter(h, v, c=colors, s=2)
        ax.set_aspect("equal")
    for k in np.unique(labels[labels >= 0]):
        cx, cy = xyz[labels == k, :2].mean(axis=0)
        axes[0].text(
            cx,
            cy,
            str(k),
            fontsize=11,
            weight="bold",
            ha="center",
            va="center",
            bbox={"facecolor": "white", "alpha": 0.7, "pad": 1, "edgecolor": "none"},
        )
    axes[0].set_title("plan (colour = plane, grey = unassigned)")
    axes[1].set_title("elevation along x")
    axes[2].set_title("elevation along y")
    fig.suptitle(title)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def save_normal_map(xyz: np.ndarray, normals: np.ndarray, out_path: str | Path, title: str = "") -> None:
    """Plan view coloured by each point's own normal: hue = compass direction the surface
    faces, saturation = pitch (white = flat, full colour >= 30 degrees). Planar facets show
    as patches of one colour; curved or twisted surfaces as gradients; ridges as thin
    bands of in-between colours. Independent of any segmentation, so it's a check on it.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import hsv_to_rgb

    az = (np.degrees(np.arctan2(normals[:, 0], normals[:, 1])) % 360) / 360
    pitch = np.degrees(np.arccos(np.clip(np.abs(normals[:, 2]), 0, 1)))
    hsv = np.column_stack([np.nan_to_num(az), np.nan_to_num(np.clip(pitch / 30, 0, 1)), np.full(len(az), 0.9)])
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(xyz[:, 0], xyz[:, 1], c=hsv_to_rgb(hsv), s=3)
    ax.set_aspect("equal")
    ax.set_title(title or "hue = facing direction, white = flat")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def write_segment_ply(xyz: np.ndarray, labels: np.ndarray, path: str | Path) -> None:
    import open3d as o3d

    pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(xyz))
    pcd.colors = o3d.utility.Vector3dVector(segment_colors(labels))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not o3d.io.write_point_cloud(str(path), pcd, write_ascii=False):
        raise OSError(f"Open3D failed to write {path}")


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
