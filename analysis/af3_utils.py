#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
af3_utils.py  (updated)

- Robustly infers seed/sample from EITHER the directory name (e.g., seed-1_sample-0)
  OR the filename (e.g., clust1_seed-1_sample-0_summary_confidences.json).
- Silently ignores top-level best summaries (without seed/sample) to avoid noise.
- Works on Windows and POSIX paths.
"""

from __future__ import annotations

import os
import re
import json
import fnmatch

from typing import Dict, List, Tuple, Any, Iterable, Optional

import numpy as np
import pandas as pd


# --------------------------- Parsing helpers ---------------------------

def _safe_load_json(path: str) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def parse_summary_file(json_path: str) -> Dict[str, Any]:
    """
    Parse one summary_confidences.json and return relevant fields.
    Returns a dict with Python-native lists for JSON-serializable fields.
    """
    data = _safe_load_json(json_path)

    return {
        "iptm": data.get("iptm", None),
        "ptm": data.get("ptm", None),
        "ranking_score": data.get("ranking_score", None),
        "fraction_disordered": data.get("fraction_disordered", None),
        "has_clash": data.get("has_clash", None),
        "chain_iptm": data.get("chain_iptm", []),
        "chain_pair_iptm": data.get("chain_pair_iptm", []),
        "chain_pair_pae_min": data.get("chain_pair_pae_min", []),
    }


_seed_sample_regex = re.compile(r"seed[-_](\d+).*sample[-_](\d+)", re.IGNORECASE)


def infer_seed_and_sample_from_dir(dir_path: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Try to infer seed/sample from the directory name (e.g., 'seed-10_sample-0').
    """
    base = os.path.basename(dir_path)
    m = _seed_sample_regex.search(base)
    if m:
        return m.group(1), m.group(2)
    return None, None


def infer_seed_and_sample_from_name(name: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Try to infer seed/sample from a filename (e.g., 'clust1_seed-10_sample-0_summary_confidences.json').
    """
    m = _seed_sample_regex.search(name)
    if m:
        return m.group(1), m.group(2)
    return None, None


def rel_dataset_label(base_dir: str, run_dir: str) -> str:
    """
    Infer the dataset/simulation label as the first path component under base_dir.
    Example:
        base_dir = /path/to/inference_output
        run_dir  = /path/to/inference_output/SimulationA/seed-1_sample-0
        -> label = "SimulationA"
    """
    rel = os.path.relpath(run_dir, base_dir)
    # First component before os.sep
    top = rel.split(os.sep)[0]
    return top


# --------------------------- Collection ---------------------------

def collect_all_runs(base_dir: str, verbose: bool = False) -> pd.DataFrame:
    """
    Walk base_dir recursively, find '*_summary_confidences.json', parse the file,
    and return a DataFrame with one row per (seed, sample) run.
    We intentionally skip summary files that do not have seed/sample in either their parent
    folder OR in the filename (to avoid top-level best summaries).
    """
    records: List[Dict[str, Any]] = []

    for root, _, files in os.walk(base_dir):
        for fname in files:
            # only look at summary files (avoid *_confidences.json)
            if not (fname == "summary_confidences.json" or fname.endswith("_summary_confidences.json")):
                continue

            json_path = os.path.join(root, fname)

            # Try dir-based inference first, then fallback to filename
            seed, sample = infer_seed_and_sample_from_dir(root)
            
            if seed is None or sample is None:
                seed, sample = infer_seed_and_sample_from_name(fname)

            if seed is None or sample is None:
                # Silently skip (these are likely dataset-level "best" summaries)
                if verbose:
                    print(f"[skip] are likely dataset-level best summaries: {json_path}")
                continue

            parsed = parse_summary_file(json_path)
            dataset = rel_dataset_label(base_dir, root)

            records.append({
                "dataset": dataset,
                "seed": seed,
                "sample": sample,
                "iptm": parsed["iptm"],
                "ptm": parsed["ptm"],
                "ranking_score": parsed["ranking_score"],
                "fraction_disordered": parsed["fraction_disordered"],
                "has_clash": parsed["has_clash"],
                # store JSON-encoded strings to keep CSV simple
                "chain_iptm": json.dumps(parsed["chain_iptm"]),
                "chain_pair_iptm": json.dumps(parsed["chain_pair_iptm"]),
                "chain_pair_pae_min": json.dumps(parsed["chain_pair_pae_min"]),
                "json_path": json_path,
                "run_folder": root,
            })

    if not records:
        return pd.DataFrame(columns=[
            "dataset","seed","sample","iptm","ptm","ranking_score",
            "fraction_disordered","has_clash",
            "chain_iptm","chain_pair_iptm","chain_pair_pae_min",
            "json_path","run_folder"
        ])

    df = pd.DataFrame.from_records(records)
    # Ensure stable column order
    cols = [
        "dataset","seed","sample","iptm","ptm","ranking_score",
        "fraction_disordered","has_clash",
        "chain_iptm","chain_pair_iptm","chain_pair_pae_min",
        "run_folder","json_path"
    ]
    df = df[[c for c in cols if c in df.columns]]
    return df


# --------------------------- Interface metrics ---------------------------

def _off_diagonal_values(mat: List[List[float]]) -> List[float]:
    """
    Return all off-diagonal entries from a square matrix 'mat'.
    If 'mat' is malformed, return [] gracefully.
    """
    try:
        n = len(mat)
        vals = []
        for i in range(n):
            row = mat[i]
            for j in range(n):
                if i == j:
                    continue
                vals.append(row[j])
        return vals
    except Exception:
        return []


def compute_iface_metrics_for_row(row: pd.Series) -> pd.Series:
    """
    From chain_pair_iptm and chain_pair_pae_min (JSON strings), compute:
        iface_pTM_max  = max off-diagonal chain-pair iptm
        iface_minPAE   = min off-diagonal chain-pair min-PAE
    Returns a Series with those two metrics.
    """
    try:
        cp_iptm = json.loads(row["chain_pair_iptm"]) if pd.notna(row["chain_pair_iptm"]) else []
        cp_pae  = json.loads(row["chain_pair_pae_min"]) if pd.notna(row["chain_pair_pae_min"]) else []
    except Exception:
        cp_iptm, cp_pae = [], []

    ptm_vals = _off_diagonal_values(cp_iptm)
    pae_vals = _off_diagonal_values(cp_pae)

    iface_pTM_max = float(np.nanmax(ptm_vals)) if len(ptm_vals) else np.nan
    iface_minPAE  = float(np.nanmin(pae_vals)) if len(pae_vals) else np.nan

    return pd.Series({"iface_pTM_max": iface_pTM_max, "iface_minPAE": iface_minPAE})


def add_interface_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Append iface metrics to a DataFrame returned by collect_all_runs.
    """
    if df.empty:
        return df.copy()
    iface = df.apply(compute_iface_metrics_for_row, axis=1)
    out = pd.concat([df.reset_index(drop=True), iface], axis=1)
    return out


# --------------------------- Plotting ---------------------------



def plot_violin_by_dataset(df: pd.DataFrame,
                           value_col: str,
                           title: str,
                           ylabel: str,
                           out_png: Optional[str] = None) -> None:
    """
    Matplotlib violin plot grouped by 'dataset' for the given value_col.
    - One chart per figure.
    - No explicit colors/styles set.
    - Orders datasets by median (descending for iPTM-like, ascending for PAE-like).
    - Shows medians inside violins.
    """
    import matplotlib.pyplot as plt

    if df.empty or value_col not in df.columns:
        raise ValueError(f"No data to plot for column '{value_col}'.")

    ascending = True if "pae" in value_col.lower() else False
    med = df.groupby("dataset")[value_col].median().sort_values(ascending=ascending)
    ordered = med.index.tolist()

    data_to_plot = [df.loc[df["dataset"] == ds, value_col].dropna().values for ds in ordered]
    med_vals = [med.loc[ds] for ds in ordered]

    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()
    parts = ax.violinplot(dataset=data_to_plot, showmeans=False, showmedians=True, showextrema=False)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(1, len(ordered) + 1))
    ax.set_xticklabels(ordered, rotation=30, ha="right")

    # Annotate median numbers slightly above the median line
    for i, m in enumerate(med_vals, start=1):
        ax.annotate(f"{m:.3f}", xy=(i, m), xytext=(0, 4),
                    textcoords="offset points", ha="center", va="bottom", fontsize=9)
    
    y_min, y_max = ax.get_ylim()
    pad = (y_max - y_min) * 0.06 if y_max > y_min else 0.1
    ax.set_ylim(y_min, y_max + pad)

    plt.tight_layout()

    if out_png:
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.show()


def plot_box_by_dataset(df: pd.DataFrame,
                        value_col: str,
                        title: str,
                        ylabel: str,
                        out_png: Optional[str] = None) -> None:
    """
    Use matplotlib to plot a simple boxplot grouped by 'dataset' for the given value_col.
    NOTE: Per instruction, we use matplotlib only, one chart per figure, and do not set colors.
    """
    import matplotlib.pyplot as plt

    if df.empty or value_col not in df.columns:
        raise ValueError(f"No data to plot for column '{value_col}'.")

    # Ensure categories sorted by median (descending for iPTM, ascending for PAE)
    ascending = True if "pae" in value_col.lower() else False
    med = df.groupby("dataset")[value_col].median().sort_values(ascending=ascending)
    ordered = med.index.tolist()

    fig = plt.figure(figsize=(10, 6))
    data_to_plot = [df.loc[df["dataset"] == ds, value_col].dropna().values for ds in ordered]

    plt.boxplot(data_to_plot, labels=ordered, showfliers=False)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    if out_png:
        fig.savefig(out_png, dpi=200, bbox_inches="tight")

    # If running in an interactive notebook, show; in scripts this is harmless.
    plt.show()


def plot_dots_by_dataset(df: pd.DataFrame,
                         value_col: str,
                         title: str,
                         ylabel: str,
                         out_png: Optional[str] = None,
                         jitter: float = 0.15,
                         s: float = 16.0) -> None:
    """
    Jittered dot plot of value_col grouped by 'dataset' (no seaborn, one chart per figure).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if df.empty or value_col not in df.columns:
        raise ValueError(f"No data to plot for column '{value_col}'.")

    ordered = sorted(df["dataset"].unique().tolist())
    x_pos = {ds: i for i, ds in enumerate(ordered, start=1)}

    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()

    for ds in ordered:
        vals = df.loc[df["dataset"] == ds, value_col].astype(float)
        x = np.full(len(vals), x_pos[ds], dtype=float) + (np.random.rand(len(vals)) - 0.5) * 2 * jitter
        ax.scatter(x, vals, s=s, alpha=0.8)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(1, len(ordered) + 1))
    ax.set_xticklabels(ordered, rotation=30, ha="right")
    plt.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.show()


# --------------------------- Summaries ---------------------------

def summarize_by_dataset(df: pd.DataFrame,
                         iptm_threshold: float = 0.7) -> pd.DataFrame:
    """
    Return a compact summary per dataset: n_runs, iptm median/IQR, fraction above threshold,
    iface_pTM_max median/IQR, iface_minPAE median/IQR.
    """
    if df.empty:
        return pd.DataFrame()

    def iqr(x: Iterable[float]) -> float:
        q75, q25 = np.nanpercentile(x, 75), np.nanpercentile(x, 25)
        return float(q75 - q25)

    g = df.groupby("dataset")
    summary = pd.DataFrame({
        "n_runs": g.size(),
        "iptm_median": g["iptm"].median(),
        "iptm_iqr": g["iptm"].apply(iqr),
        "iptm_frac_ge_thresh": g["iptm"].apply(lambda s: float(np.mean(s >= iptm_threshold))),
        "iface_pTM_max_median": g["iface_pTM_max"].median(),
        "iface_pTM_max_iqr": g["iface_pTM_max"].apply(iqr),
        "iface_minPAE_median": g["iface_minPAE"].median(),
        "iface_minPAE_iqr": g["iface_minPAE"].apply(iqr),
    }).reset_index()

    return summary