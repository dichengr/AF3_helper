#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
af3_utils.py (Updated)

Utilities for collecting, processing, and plotting AlphaFold3 results.
- Integrates file finding (PAE JSON, CIF) into the initial run collection.
- Uses flexible file patterns inspired by ipsae.py for robustness.
"""

from __future__ import annotations

import os
import re
import json
import fnmatch
import subprocess
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import pandas as pd

# --------------------------- Parsing & File Finding ---------------------------

def _safe_load_json(path: str) -> dict:
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"[-] Warning: Could not read or parse JSON {path}: {e}")
        return {}

_seed_sample_regex = re.compile(r"seed[-_](\d+).*sample[-_](\d+)", re.IGNORECASE)

def _infer_seed_and_sample(path: str, filename: str) -> Tuple[Optional[str], Optional[str]]:
    """Try to infer seed/sample from directory path first, then filename."""
    dir_match = _seed_sample_regex.search(os.path.basename(path))
    if dir_match:
        return dir_match.group(1), dir_match.group(2)
    file_match = _seed_sample_regex.search(filename)
    if file_match:
        return file_match.group(1), file_match.group(2)
    return None, None

def rel_dataset_label(base_dir: str, run_dir: str) -> str:
    """Infer the dataset label from the first path component under base_dir."""
    rel = os.path.relpath(run_dir, base_dir)
    return rel.split(os.sep)[0]

    # More flexible patterns inspired by ipsae.py

_AF3_JSON_PATTERNS = ["*confidences.json"]
_AF3_CIF_PATTERNS = ["*.cif"]

def _find_first_match(directory: str, patterns: list[str], exclude: str = "") -> str | None:
    """Find the first file in a directory matching a list of patterns."""
    for pat in patterns:
        for fn in sorted(fnmatch.filter(os.listdir(directory), pat)):
            if exclude and exclude in fn:
                continue
            return os.path.join(directory, fn)
    return None

def find_ipsae_inputs(run_folder: str) -> Tuple[Optional[str], Optional[str]]:
    """Return (pae_json_path, cif_path) for a run folder using flexible patterns."""
    # For JSON, we must exclude 'summary_confidences.json'
    json_path = _find_first_match(run_folder, _AF3_JSON_PATTERNS, exclude="summary")
    cif_path = _find_first_match(run_folder, _AF3_CIF_PATTERNS)
    return json_path, cif_path

# --------------------------- Collection (Integrated) ---------------------------

def collect_all_runs(base_dir: str) -> pd.DataFrame:
    """
    Walk base_dir, find valid runs, parse summaries, and find input files for ipSAE.
    """
    records: List[Dict[str, Any]] = []
    for root, _, files in os.walk(base_dir):
        # Find summary files to identify valid runs
        summary_files = [f for f in files if "summary_confidences.json" in f]
        if not summary_files:
            continue

        for fname in summary_files:
            seed, sample = _infer_seed_and_sample(root, fname)
            if seed is None or sample is None:
                continue  # Skip top-level or non-run-specific summaries

            # A valid run was found, now parse its data and find its input files
            summary_path = os.path.join(root, fname)
            parsed_summary = _safe_load_json(summary_path)
            pae_json_path, cif_path = find_ipsae_inputs(root)

            records.append({
                "dataset": rel_dataset_label(base_dir, root),
                "seed": seed,
                "sample": sample,
                "iptm": parsed_summary.get("iptm"),
                "ptm": parsed_summary.get("ptm"),
                "ranking_score": parsed_summary.get("ranking_score"),
                "chain_pair_iptm": json.dumps(parsed_summary.get("chain_pair_iptm", [])),
                "chain_pair_pae_min": json.dumps(parsed_summary.get("chain_pair_pae_min", [])),
                "run_folder": root,
                "summary_path": summary_path,
                "cif_path": cif_path,
                "pae_json_path": pae_json_path,
            })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame.from_records(records)
    # Ensure stable column order
    cols = [
        "dataset", "seed", "sample", "iptm", "ptm", "ranking_score",
        "chain_pair_iptm", "chain_pair_pae_min", "run_folder", "cif_path",
        "pae_json_path", "summary_path"
    ]
    return df[cols]

# --------------------------- Interface & Plotting --------------------

def _off_diagonal_values(mat: List[List[float]]) -> List[float]:
    if not isinstance(mat, list) or not all(isinstance(row, list) for row in mat):
        return []
    n = len(mat)

    return [mat[i][j] for i in range(n) for j in range(n) if i != j]

def compute_iface_metrics_for_row(row: pd.Series) -> pd.Series:
    try:
        cp_iptm = json.loads(row["chain_pair_iptm"]) if pd.notna(row["chain_pair_iptm"]) else []
        cp_pae = json.loads(row["chain_pair_pae_min"]) if pd.notna(row["chain_pair_pae_min"]) else []
    except (json.JSONDecodeError, TypeError):
        cp_iptm, cp_pae = [], []

    ptm_vals = _off_diagonal_values(cp_iptm)
    pae_vals = _off_diagonal_values(cp_pae)
    iface_pTM_max = float(np.nanmax(ptm_vals)) if ptm_vals else np.nan
    iface_minPAE = float(np.nanmin(pae_vals)) if pae_vals else np.nan

    return pd.Series({"iface_pTM_max": iface_pTM_max, "iface_minPAE": iface_minPAE})

def add_interface_metrics(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    
    iface_metrics = df.apply(compute_iface_metrics_for_row, axis=1)

    return pd.concat([df.reset_index(drop=True), iface_metrics], axis=1)

def plot_by_dataset(df: pd.DataFrame, value_col: str, title: str, ylabel: str,
                    plot_type: str = 'dots', out_png: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt

    if df.empty or value_col not in df.columns or df[value_col].isnull().all():
        print(f"[-] Warning: No data to plot for column '{value_col}'. Skipping plot.")
        return
    ascending = "pae" in value_col.lower()
    med = df.groupby("dataset")[value_col].median().sort_values(ascending=ascending)
    ordered = med.index.tolist()
    data_to_plot = [df.loc[df["dataset"] == ds, value_col].dropna().values for ds in ordered]
    fig, ax = plt.subplots(figsize=(10, 6))

    if plot_type == 'box':
        ax.boxplot(data_to_plot, showfliers=False)
    elif plot_type == 'violin':
        ax.violinplot(dataset=data_to_plot, showmedians=True)
    elif plot_type == 'dots':
        for i, ds_data in enumerate(data_to_plot, start=1):
            x = np.random.normal(i, 0.08, size=len(ds_data))
            ax.scatter(x, ds_data, s=16, alpha=0.6)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(1, len(ordered) + 1))
    ax.set_xticklabels(ordered, rotation=30, ha="right")
    plt.tight_layout()
    
    if out_png:
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.show()
    plt.close(fig)

# --------------------------- ipSAE Integration (Simplified) --------------------

def run_ipsae_on_run(run_folder: str, pae_json_path: str, cif_path: str, ipsae_py: str,
                       pae_cutoff: int, dist_cutoff: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Execute ipsae.py on a single run using pre-found file paths."""
    if not pae_json_path or not cif_path:
        print(f"[-] Warning: Missing PAE JSON or CIF path for {run_folder}. Skipping IPSAE.")
        return pd.DataFrame(), pd.DataFrame()

    cif_stem = os.path.basename(cif_path).replace(".cif", "")
    pae_str, dist_str = f"{pae_cutoff:02d}", f"{dist_cutoff:02d}"
    base_name = f"{cif_stem}_{pae_str}_{dist_str}"
    chain_output_path = os.path.join(run_folder, f"{base_name}.txt")
    res_output_path = os.path.join(run_folder, f"{base_name}_byres.txt")

    cmd = ["python", ipsae_py, pae_json_path, cif_path, str(pae_cutoff), str(dist_cutoff)]
    proc = subprocess.run(cmd, cwd=run_folder, capture_output=True, text=True, check=False)

    if proc.returncode != 0:
        print(f"[!] Error running IPSAE on {run_folder}:\n{proc.stderr}")
        return pd.DataFrame(), pd.DataFrame()

    df_chain = pd.read_csv(chain_output_path, sep=r"\s+", on_bad_lines='skip').dropna(how='all') if os.path.exists(chain_output_path) else pd.DataFrame()
    df_res = pd.read_csv(res_output_path, sep=r"\s+", on_bad_lines='skip').dropna(how='all') if os.path.exists(res_output_path) else pd.DataFrame()
    return df_chain, df_res

def collect_ipsae_for_runs(df_runs: pd.DataFrame, ipsae_py: str, pae_cutoff: int,
                           dist_cutoff: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """For each run in the DataFrame, run IPSAE and concatenate the results."""
    if df_runs.empty:
        return pd.DataFrame(), pd.DataFrame()

    all_chains, all_res = [], []
    for _, row in df_runs.iterrows():
        cdf, rdf = run_ipsae_on_run(
            run_folder=row["run_folder"],
            pae_json_path=row["pae_json_path"],
            cif_path=row["cif_path"],
            ipsae_py=ipsae_py,
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff
        )
        if not cdf.empty:
            for col in ["dataset", "seed", "sample"]: cdf[col] = row[col]
            all_chains.append(cdf)
        if not rdf.empty:
            for col in ["dataset", "seed", "sample"]: rdf[col] = row[col]
            all_res.append(rdf)

    return pd.concat(all_chains, ignore_index=True) if all_chains else pd.DataFrame(), \
           pd.concat(all_res, ignore_index=True) if all_res else pd.DataFrame()