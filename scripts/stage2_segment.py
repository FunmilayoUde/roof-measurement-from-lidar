"""Stage 2: split one house's roof points into planes and measure each plane.

    uv run python scripts/stage2_segment.py --house-id delft_01 --method rg
    uv run python scripts/stage2_segment.py --house-id delft_01 --method ransac --threshold 0.08

Reads data/interim/<house-id>/{roof.ply, stats.json} from stage 1. Writes
data/interim/<house-id>/stage2_<method>/{segments.ply, planes.json, preview.png}
and prints the row to paste into RESULTS.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import open3d as o3d

from roofrecon.report import summarize
from roofrecon.segment import segment_ransac, segment_region_growing
from roofrecon.normals import estimate_normals
from roofrecon.viz import save_normal_map, save_segment_snapshot, write_segment_ply


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--house-id", required=True)
    ap.add_argument("--method", choices=["ransac", "rg"], required=True)
    ap.add_argument(
        "--threshold",
        type=float,
        default=0.12,
        help="Max point-to-plane distance, m (chosen by scripts/stage2_sweep.py)",
    )
    ap.add_argument(
        "--min-area",
        type=float,
        default=1.0,
        help="Smallest plane to keep, m² of plan area; converted to points with the house's density",
    )
    ap.add_argument("--connect-radius", type=float, default=0.3, help="Connectivity / grow radius, m")
    ap.add_argument("--normal-radius", type=float, default=0.5, help="rg only: neighbourhood for normals, m")
    ap.add_argument("--angle", type=float, default=8.0, help="rg only: max normal deviation from region, degrees")
    ap.add_argument("--out-suffix", default="", help="Appended to the output folder name, for sweeps")
    ap.add_argument("--root", type=Path, default=Path("data/interim"))
    args = ap.parse_args()

    house = args.root / args.house_id
    xyz = np.asarray(o3d.io.read_point_cloud(str(house / "roof.ply")).points)
    stats = json.loads((house / "stats.json").read_text())
    density = stats["density_inside_walls"]["nominal_pts_per_m2"]
    min_points = max(3, int(round(args.min_area * density)))

    if args.method == "ransac":
        seg = segment_ransac(xyz, threshold=args.threshold, min_points=min_points, connect_radius=args.connect_radius)
    else:
        seg = segment_region_growing(
            xyz,
            normal_radius=args.normal_radius,
            angle_deg=args.angle,
            threshold=args.threshold,
            grow_radius=args.connect_radius,
            min_points=min_points,
        )

    summary = {
        "house_id": args.house_id,
        "method": args.method,
        "params": {k: v for k, v in vars(args).items() if k not in ("root", "house_id", "method", "out_suffix")}
        | {"min_points": min_points, "density_pts_per_m2": density},
    } | summarize(xyz, seg, density)
    planes = summary["per_plane"]
    roofs = [d for d in planes if d["kind"] != "wall"]

    out = house / f"stage2_{args.method}{args.out_suffix}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "planes.json").write_text(json.dumps(summary, indent=2))
    write_segment_ply(xyz, seg.labels, out / "segments.ply")
    normals, _, _ = estimate_normals(xyz, args.normal_radius)
    save_normal_map(xyz, normals, out / "normal_map.png", title=f"{args.house_id}: hue = facing, white = flat")
    save_segment_snapshot(
        xyz,
        seg.labels,
        out / "preview.png",
        title=f"{args.house_id} {args.method}: {len(planes)} planes, coverage {seg.coverage:.1%}, "
        f"threshold {args.threshold} m",
    )

    print(f"{args.house_id} [{args.method}] min_points={min_points} (min_area {args.min_area} m² × {density:.1f}/m²)")
    print(
        f"  {'id':>3} {'pts':>5} {'kind':>6} {'pitch°':>7} {'az°':>6} {'rms m':>7} "
        f"{'true m²':>8} {'plan m²':>8} {'fill':>5}"
    )
    for d in planes:
        az = "-" if d["azimuth_deg"] is None else f"{d['azimuth_deg']:.0f}"
        print(
            f"  {d['id']:>3} {d['points']:>5} {d['kind']:>6} {d['pitch_deg']:>7.1f} {az:>6} {d['rms_m']:>7.3f} "
            f"{d['true_area_m2']:>8.1f} {d['plan_area_m2']:>8.1f} {d['fill']:>5.2f}"
        )
    print("\nRESULTS.md row:")
    print(
        f"| {args.house_id} | {args.method} | {args.threshold} | {min_points} | {len(planes)} | {len(roofs)} | "
        f"{summary['coverage']:.3f} | {summary['median_rms_m']} | {summary['roof_true_area_m2']} |"
    )


if __name__ == "__main__":
    main()
