"""
af3_one_shot.py

One-shot pipeline:
    Input  : a BASE inference-output folder that contains multiple simulation subfolders
    Output : (1) all_runs_summary.csv (raw per-run metrics)
             (2) all_interface_metrics.csv (with iface_pTM_max & iface_minPAE)
             (3) boxplots comparing datasets on iptm and iface_minPAE

Usage:
    python af3_one_shot.py --base_dir /path/to/inference_output --out_dir  /path/to/save_results   --iptm_thresh 0.7

Notes:
  * We skip any summary_confidences.json that are not inside a folder named
    'seed-*_sample-*' to avoid collecting top-level "best" summaries without seed/sample.
  * Boxplots are saved into out_dir as PNGs.
"""

from __future__ import annotations

import os
import argparse
import pandas as pd

from af3_utils import (
    collect_all_runs,
    add_interface_metrics,
    plot_box_by_dataset,
    plot_violin_by_dataset,
    plot_dots_by_dataset,
    summarize_by_dataset,
)


def parse_args():
    p = argparse.ArgumentParser(description="Collect AF3 metrics across many simulations and plot comparisons.")
    p.add_argument("--base_dir", required=True, help="Root folder containing simulation subfolders.")
    p.add_argument("--out_dir", default=None, help="Where to write CSVs/plots (default: BASE_DIR).")
    p.add_argument("--iptm_thresh", type=float, default=0, help="Threshold for robust-run fraction summary (default: 0.7).")
    # NEW plot style switch (default: dots for early-stage)
    p.add_argument("--plot", choices=["dots","box","violin"], default="dots",
                   help="Style for per-dataset plots of key metrics (default: dots).")
    
    return p.parse_args()


def main():
    args = parse_args()
    base_dir = os.path.abspath(args.base_dir)
    out_dir  = os.path.abspath(args.out_dir or base_dir)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[+] Scanning under: {base_dir}")
    df_raw = collect_all_runs(base_dir, verbose=True)
    if df_raw.empty:
        raise SystemExit("No runs found (no '*summary_confidences.json' under seed-*_sample-* folders).")

    raw_csv = os.path.join(out_dir, "all_runs_summary.csv")
    df_raw.to_csv(raw_csv, index=False)
    print(f"[+] Wrote raw runs CSV: {raw_csv}  (n={len(df_raw)})")

    print("[+] Computing interface metrics...")
    df_iface = add_interface_metrics(df_raw)
    iface_csv = os.path.join(out_dir, "all_interface_metrics.csv")
    df_iface.to_csv(iface_csv, index=False)
    print(f"[+] Wrote runs+iface CSV: {iface_csv}")

    # Summaries per dataset
    df_sum = summarize_by_dataset(df_iface, iptm_threshold=args.iptm_thresh)
    sum_csv = os.path.join(out_dir, "dataset_summaries.csv")
    df_sum.to_csv(sum_csv, index=False)
    print(f"[+] Wrote dataset summaries: {sum_csv}")


    # If you still want the legacy plots, you can conditionally call:
    if args.plot == "dots":
        # iPTM dots
        iptm_png = os.path.join(out_dir, "dots_iptm_by_dataset.png")
        plot_dots_by_dataset(df_iface, value_col="iptm",
                             title="iPTM by dataset", ylabel="iPTM",
                             out_png=iptm_png)
        print(f"[+] Saved: {iptm_png}")

        pae_png = os.path.join(out_dir, "dots_iface_minPAE_by_dataset.png")
        plot_dots_by_dataset(df_iface, value_col="iface_minPAE", title="Interface min PAE by dataset", ylabel="iface_minPAE (Å)", out_png=pae_png)
        print(f"[+] Saved: {pae_png}")
        
    elif args.plot == "box":
        iptm_png = os.path.join(out_dir, "box_iptm_by_dataset.png")
        plot_box_by_dataset(df_iface, value_col="iptm", title="iPTM by dataset", ylabel="iPTM", out_png=iptm_png)
        print(f"[+] Saved: {iptm_png}")

        pae_png = os.path.join(out_dir, "box_iface_minPAE_by_dataset.png")
        plot_box_by_dataset(df_iface, value_col="iface_minPAE", title="Interface min PAE by dataset", ylabel="iface_minPAE (Å)", out_png=pae_png)
        print(f"[+] Saved: {pae_png}")

    elif args.plot == "violin":
        iptm_png = os.path.join(out_dir, "violin_iptm_by_dataset.png")
        plot_violin_by_dataset(df_iface, value_col="iptm", title="iPTM by dataset", ylabel="iPTM", out_png=iptm_png)
        print(f"[+] Saved: {iptm_png}")

        pae_png = os.path.join(out_dir, "violin_iface_minPAE_by_dataset.png")
        plot_violin_by_dataset(df_iface, value_col="iface_minPAE", title="Interface min PAE by dataset", ylabel="iface_minPAE (Å)", out_png=pae_png)
        print(f"[+] Saved: {pae_png}")

    print("[✓] Done.")


if __name__ == "__main__":
    main()
