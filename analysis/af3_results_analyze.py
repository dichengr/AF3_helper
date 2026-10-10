"""
One-shot analysis pipeline for AlphaFold3 inference outputs.
"""

from __future__ import annotations

import os
import argparse
import pandas as pd

from af3_utils import (
    collect_all_runs,
    add_interface_metrics,
    collect_ipsae_for_runs,
    plot_by_dataset
)

def parse_args():
    p = argparse.ArgumentParser(description="Collect and analyze AF3 metrics across multiple simulation runs.")
    p.add_argument("--base_dir", required=True, help="Root folder containing simulation subfolders.")
    p.add_argument("--out_dir", required=True, help="Directory to save output CSVs and plots.")
    p.add_argument("--ipsae_py", default=r'E:\University\RA\Af3_helper\ipsae.py', help="Path to ipsae.py to enable advanced metric calculation.")
    p.add_argument("--ipsae_pae", type=int, default=10, help="ipSAE PAE cutoff (default: 10).")
    p.add_argument("--ipsae_dist", type=int, default=10, help="ipSAE distance cutoff (default: 10).")
    p.add_argument("--plot", choices=["dots", "box", "violin"], default="dots", help="Plot style for metrics (default: dots).")
    return p.parse_args()

def main():
    args = parse_args()
    base_dir = os.path.abspath(args.base_dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[1/4] Scanning for runs and input files under: {base_dir}")
    df_raw = collect_all_runs(base_dir)
    if df_raw.empty:
        raise SystemExit("[!] No valid runs found. Ensure folders are named '...seed-X_sample-Y...'.")
    print(f"    Found {len(df_raw)} runs across {df_raw['dataset'].nunique()} datasets.")

    print("[2/4] Computing interface metrics (iface_pTM_max, iface_minPAE)...")
    df_iface = add_interface_metrics(df_raw)
    iface_csv = os.path.join(out_dir, "all_interface_metrics.csv")
    df_iface.to_csv(iface_csv, index=False)
    print(f"    -> Saved main results to: {iface_csv}")

    # Plot standard metrics
    plot_by_dataset(
        df_iface, "iptm", "iPTM by Dataset", "iPTM Score",
        plot_type=args.plot, out_png=os.path.join(out_dir, f"{args.plot}_iptm_by_dataset.png")
    )
    plot_by_dataset(
        df_iface, "iface_minPAE", "Interface min PAE by Dataset", "min PAE (Å)",
        plot_type=args.plot, out_png=os.path.join(out_dir, f"{args.plot}_iface_minPAE_by_dataset.png")
    )

    if args.ipsae_py:
        print(f"[3/4] Running `ipsae.py` on {len(df_iface)} runs...")
        chain_df, res_df = collect_ipsae_for_runs(
            df_runs=df_iface,
            ipsae_py=args.ipsae_py,
            pae_cutoff=args.ipsae_pae,
            dist_cutoff=args.ipsae_dist
        )
        if not chain_df.empty:
            ipsae_chain_csv = os.path.join(out_dir, "ipsae_chainpair_scores.csv")
            chain_df.to_csv(ipsae_chain_csv, index=False)
            print(f"    -> Saved ipSAE chain-pair scores to: {ipsae_chain_csv}")

            if not res_df.empty:
                ipsae_res_csv = os.path.join(out_dir, "ipsae_residue_scores.csv")
                res_df.to_csv(ipsae_res_csv, index=False)
                print(f"    -> Saved ipSAE residue scores to: {ipsae_res_csv}")

            print("[4/4] Generating plots for advanced metrics...")
            plot_df = chain_df.loc[chain_df["Type"].str.lower() == "max"].copy() if "Type" in chain_df.columns else chain_df.copy()
            metrics_to_plot = {"ipSAE": "ipSAE Score", "pDockQ": "pDockQ Score", "pDockQ2": "pDockQ2 Score", "LIS": "LIS Score"}

            for col, label in metrics_to_plot.items():
                if col in plot_df.columns:
                    plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")
                    plot_by_dataset(
                        plot_df, col, f"{col} (Type=max) by Dataset", label,
                        plot_type=args.plot, out_png=os.path.join(out_dir, f"{args.plot}_{col}_by_dataset.png")
                    )
        else:
            print("    -> No results generated from `ipsae.py`.")
    else:
        print("[3/4] Skipping `ipsae.py` run (path not provided).")
        print("[4/4] Skipping advanced metric plotting.")

    print("\n[✓] Analysis complete.")

if __name__ == "__main__":
    main()