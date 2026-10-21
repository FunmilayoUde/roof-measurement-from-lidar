"""Stage 2 parameter sweeps: vary one setting at a time, everything else at its default.

    uv run python scripts/stage2_sweep.py

Writes data/interim/sweeps/stage2_sweep.csv (one row per house x method x setting value)
and docs/figures/sweep_<setting>.png. Reference checks:
  delft_02 must give exactly 2 planes, both flat;
  delft_03 must give exactly 5 planes, with its two gable halves within 5% in area.
delft_01 (curved top) and delft_05 (ambiguous) are reported but don't decide settings.
"""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
import open3d as o3d

from roofrecon.report import summarize
from roofrecon.segment import segment_ransac, segment_region_growing

ROOT = Path("data/interim")
HOUSES = ["delft_02", "delft_03", "delft_01", "delft_05"]
REFERENCE = {"delft_02", "delft_03"}
DEFAULTS = {"threshold": 0.10, "min_area": 1.0, "normal_radius": 0.5, "angle": 8.0}
SWEEPS = {
    "threshold": ([0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25], ["ransac", "rg"]),
    "min_area": ([0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0], ["ransac", "rg"]),
    "normal_radius": ([0.3, 0.4, 0.5, 0.6, 0.75, 1.0], ["rg"]),
    "angle": ([4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0], ["rg"]),
}
# Picked from the joint check: centre of the region where every seed passes both references.
CHOSEN = {"threshold": 0.12, "min_area": 1.0}
LABELS = {"threshold": "distance threshold (m)", "min_area": "minimum plane area (m²)",
          "normal_radius": "normal radius (m)", "angle": "angle tolerance (°)"}


def run(xyz: np.ndarray, density: float, method: str, p: dict):
    min_points = max(3, int(round(p["min_area"] * density)))
    if method == "ransac":
        return segment_ransac(xyz, threshold=p["threshold"], min_points=min_points)
    return segment_region_growing(xyz, normal_radius=p["normal_radius"], angle_deg=p["angle"],
                                  threshold=p["threshold"], min_points=min_points)


def reference_ok(house: str, s: dict) -> bool | None:
    planes = s["per_plane"]
    if house == "delft_02":
        return s["planes"] == 2 and all(d["kind"] == "flat" for d in planes)
    if house == "delft_03":
        sloped = sorted((d for d in planes if d["kind"] == "sloped"), key=lambda d: -d["points"])[:2]
        if s["planes"] != 5 or len(sloped) < 2:
            return False
        a, b = sloped[0]["true_area_m2"], sloped[1]["true_area_m2"]
        return abs(a - b) / max(a, b) <= 0.05
    return None


def main() -> None:
    data = {}
    for h in HOUSES:
        xyz = np.asarray(o3d.io.read_point_cloud(str(ROOT / h / "roof.ply")).points)
        density = json.loads((ROOT / h / "stats.json").read_text())["density_inside_walls"]["nominal_pts_per_m2"]
        data[h] = (xyz, density)

    rows = []
    for setting, (values, methods) in SWEEPS.items():
        for v in values:
            params = DEFAULTS | {setting: v}
            for method in methods:
                for h in HOUSES:
                    xyz, density = data[h]
                    t0 = time.time()
                    s = summarize(xyz, run(xyz, density, method, params), density)
                    rows.append({
                        "setting": setting, "value": v, "method": method, "house": h,
                        "planes": s["planes"], "coverage": s["coverage"], "median_rms_m": s["median_rms_m"],
                        "roof_true_area_m2": s["roof_true_area_m2"], "reference_ok": reference_ok(h, s),
                        "seconds": round(time.time() - t0, 2),
                    })
                    print(f"{setting}={v:<5} {method:6} {h}: planes {s['planes']:>2} cov {s['coverage']:.3f} "
                          f"rms {s['median_rms_m']} area {s['roof_true_area_m2']} ref {rows[-1]['reference_ok']}")

    out = ROOT / "sweeps"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "stage2_sweep.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for setting in SWEEPS:
        plot_sweep([r for r in rows if r["setting"] == setting], setting, Path("docs/figures") / f"sweep_{setting}.png")
    print(f"\n{len(rows)} runs -> {out / 'stage2_sweep.csv'}")
    joint_ransac_check(data, out)


def joint_ransac_check(data: dict, out: Path, seeds: range = range(10)) -> None:
    """One-at-a-time sweeps use one RANSAC seed and vary one setting. This checks the
    threshold x min_area combinations over several seeds: the fraction of seeds for which
    both reference houses pass, then plane-count stability per house at the chosen pair."""
    thresholds, min_areas = [0.08, 0.10, 0.12, 0.15, 0.20], [0.75, 1.0, 1.5, 2.0]
    grid = []
    print("\nRANSAC: fraction of seeds passing both references")
    print("threshold \\ min_area " + "".join(f"{m:>7}" for m in min_areas))
    for t in thresholds:
        line = []
        for m in min_areas:
            passed = 0
            for seed in seeds:
                ok = True
                for h in sorted(REFERENCE):
                    xyz, d = data[h]
                    seg = segment_ransac(xyz, threshold=t, min_points=max(3, int(round(m * d))), seed=seed)
                    ok &= bool(reference_ok(h, summarize(xyz, seg, d)))
                passed += ok
            grid.append({"threshold": t, "min_area": m, "seeds": len(seeds), "pass_fraction": passed / len(seeds)})
            line.append(passed / len(seeds))
        print(f"{t:<21} " + "".join(f"{v:>7.0%}" for v in line))
    with (out / "stage2_joint_ransac.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(grid[0]))
        w.writeheader()
        w.writerows(grid)

    t, m = CHOSEN["threshold"], CHOSEN["min_area"]
    print(f"\nRANSAC seed stability at threshold {t}, min_area {m}:")
    for h in HOUSES:
        xyz, d = data[h]
        res = [summarize(xyz, segment_ransac(xyz, threshold=t, min_points=max(3, int(round(m * d))), seed=s), d) for s in seeds]
        counts = {v: [r["planes"] for r in res].count(v) for v in sorted({r["planes"] for r in res})}
        cov = [r["coverage"] for r in res]
        area = [r["roof_true_area_m2"] for r in res]
        print(f"  {h}: planes {counts}, coverage {min(cov):.3f}-{max(cov):.3f}, area {min(area):.1f}-{max(area):.1f} m²")


# Validated categorical slots 1-4 (reference palette, light mode), fixed per house.
HOUSE_COLORS = {"delft_02": "#2a78d6", "delft_03": "#eb6834", "delft_01": "#1baf7a", "delft_05": "#eda100"}
METHOD_STYLE = {"ransac": "-", "rg": "--"}
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"


def label_line_ends(ax, rows: list[dict], method: str, key: str) -> None:
    """House names at the right end of each line, nudged apart so they don't overlap."""
    ends = []
    for h in HOUSES:
        last = max((r for r in rows if r["house"] == h and r["method"] == method), key=lambda r: r["value"])
        ends.append([last[key], h, last["value"]])
    ends.sort()
    lo, hi = ax.get_ylim()
    gap = 0.07 * (hi - lo)
    for i in range(1, len(ends)):
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
    for y, h, x in ends:
        ax.annotate(h, (x, y), xytext=(6, 0), textcoords="offset points", color=INK_2, fontsize=8, va="center",
                    annotation_clip=False)


def plot_sweep(rows: list[dict], setting: str, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    metrics = [("planes", "planes found"), ("coverage", "coverage (fraction assigned)"),
               ("median_rms_m", "median plane rms (m)"), ("roof_true_area_m2", "roof true area (m²)")]
    methods = sorted({r["method"] for r in rows}, key=["ransac", "rg"].index)
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), facecolor=SURFACE)
    for ax, (key, title) in zip(axes, metrics):
        ax.set_facecolor(SURFACE)
        for h in HOUSES:
            for m in methods:
                pts = sorted((r["value"], r[key]) for r in rows if r["house"] == h and r["method"] == m)
                xs, ys = zip(*pts)
                ax.plot(xs, ys, METHOD_STYLE[m], color=HOUSE_COLORS[h], lw=2, marker="o", ms=4)
                # Reference failures: hollow ring on the point so they read without colour.
                if key == "planes" and h in REFERENCE:
                    bad = [(r["value"], r[key]) for r in rows if r["house"] == h and r["method"] == m and r["reference_ok"] is False]
                    if bad:
                        ax.scatter(*zip(*bad), s=90, facecolors="none", edgecolors=INK, linewidths=1.2, zorder=5)
        if key == "planes":
            label_line_ends(ax, rows, methods[0], key)
        if DEFAULTS[setting] is not None:
            ax.axvline(DEFAULTS[setting], color=INK_2, lw=1, ls=":", zorder=0)
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        ax.set_xlabel(LABELS[setting], color=INK_2, fontsize=9)
        ax.tick_params(colors=INK_2, labelsize=8)
        ax.grid(True, color=GRID, lw=0.8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
    handles = [Line2D([], [], color=HOUSE_COLORS[h], lw=2, label=h + (" (reference)" if h in REFERENCE else "")) for h in HOUSES]
    handles += [Line2D([], [], color=INK_2, lw=2, ls=METHOD_STYLE[m], label={"ransac": "RANSAC", "rg": "region growing"}[m]) for m in methods]
    handles += [Line2D([], [], ls="none", marker="o", ms=9, mfc="none", mec=INK, label="reference check failed"),
                Line2D([], [], color=INK_2, lw=1, ls=":", label="current default")]
    fig.suptitle(f"Stage 2 sweep: {LABELS[setting]} (other settings at defaults)", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8.5, labelcolor=INK)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
