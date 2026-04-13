#!/usr/bin/env python3
"""
Differential Enrichment Analysis for RewindPipeline
=====================================================
Compares clone fold-changes between Treated and Control conditions relative
to a T0 baseline to identify clones significantly enriched or depleted by
drug selection, while controlling for random drift.

Statistical framework
---------------------
1. For every clone observed at T0, compute:
       FC_control = (Control_size + pseudocount) / (T0_size + pseudocount)
       FC_treated = (Treated_size + pseudocount) / (T0_size + pseudocount)
2. Compute Log2 fold-change ratio = log2(FC_treated / FC_control).
3. Perform Fisher's exact test on [[T0_ctrl, T1_ctrl], [T0_treat, T1_treat]].
4. Apply Benjamini-Hochberg FDR correction.
5. Classify clones as Enriched / Depleted / Neutral.

Outputs
-------
- differential_enrichment_results.txt   Main results table
- enriched_clones.txt                    Enriched clones only
- depleted_clones.txt                    Depleted clones only
- neutral_clones.txt                     Neutral clones only
- summary_statistics.txt                 Overall statistics
- plots/                                 Publication-quality visualisations

Usage
-----
    python differential_enrichment_analysis.py \\
        --t0-clones T0_cloneID_cloneBarcode.tsv \\
        --control-clones Control_cloneID_cloneBarcode.tsv \\
        --treated-clones Treated_cloneID_cloneBarcode.tsv \\
        --control-matches control_matching/best_matches.txt \\
        --treated-matches treated_matching/best_matches.txt \\
        --output-dir enrichment_results
"""

import argparse
import os
import sys
import warnings
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RESULT_COLUMNS = [
    "Clone_ID", "T0_size", "Control_size", "Treated_size",
    "FC_control", "FC_treated", "Log2_FC_ratio",
    "p_value", "q_value", "Status",
]


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

def parse_clone_barcode_file(path: str) -> dict:
    """
    Parse a cloneID_cloneBarcode.tsv file.

    Returns
    -------
    dict
        {clone_id (int): set of lineage barcodes}
    """
    df = pd.read_csv(path, sep="\t")
    clone_to_barcodes: dict[int, set] = {}
    for _, row in df.iterrows():
        cid = int(row["cloneID"])
        barcodes = set(str(row["clone_barcode"]).split("_"))
        clone_to_barcodes[cid] = barcodes
    return clone_to_barcodes


def load_clone_sizes_from_filtered(path: str) -> dict:
    """
    Load total UMI counts per clone from a rewind_filtered.tsv file.

    The file has columns: cell_barcode, corrected_lineage_barcode, umi_count, ...
    Clone sizes are obtained by summing UMI counts across all cells and barcodes
    belonging to each clone (approximated from the clone barcode file).

    Returns
    -------
    dict
        {lineage_barcode: total_umi_count}
    """
    df = pd.read_csv(path, sep="\t")
    umi_col = "umi_count" if "umi_count" in df.columns else df.columns[2]
    bc_col = "corrected_lineage_barcode" if "corrected_lineage_barcode" in df.columns else df.columns[1]
    weights = df.groupby(bc_col)[umi_col].sum().to_dict()
    return weights


def compute_clone_sizes(clone_barcodes: dict, barcode_umis: dict) -> dict:
    """
    Compute total UMI count per clone by summing the UMI counts of its
    constituent lineage barcodes.

    Parameters
    ----------
    clone_barcodes : dict
        {clone_id: set of lineage barcodes}
    barcode_umis : dict
        {lineage_barcode: total_umi_count}

    Returns
    -------
    dict
        {clone_id: total_umi_count}
    """
    sizes = {}
    for cid, bcs in clone_barcodes.items():
        sizes[cid] = sum(barcode_umis.get(bc, 0) for bc in bcs)
    return sizes


def load_best_matches(path: str) -> dict:
    """
    Load best_matches.txt from match_timepoints.py.

    Returns
    -------
    dict
        {T1_clone_id: T0_clone_id}
    """
    df = pd.read_csv(path, sep="\t")
    mapping = {}
    for _, row in df.iterrows():
        t1 = int(row["T1_clone_id"])
        t0 = int(row["T0_clone_id"])
        mapping[t1] = t0
    return mapping


def build_t0_to_t1_size_map(
    t1_clone_sizes: dict,
    t1_to_t0_map: dict,
) -> dict:
    """
    Build a mapping from T0 clone IDs to T1 clone sizes using the
    timepoint matching results.

    If multiple T1 clones map to the same T0 clone, their sizes are summed.

    Parameters
    ----------
    t1_clone_sizes : dict
        {t1_clone_id: size}
    t1_to_t0_map : dict
        {t1_clone_id: t0_clone_id}

    Returns
    -------
    dict
        {t0_clone_id: aggregated_t1_size}
    """
    t0_sizes: dict[int, float] = defaultdict(float)
    for t1_cid, t0_cid in t1_to_t0_map.items():
        t0_sizes[t0_cid] += t1_clone_sizes.get(t1_cid, 0)
    return dict(t0_sizes)


# ---------------------------------------------------------------------------
# Statistical analysis
# ---------------------------------------------------------------------------

def run_enrichment_analysis(
    t0_sizes: dict,
    control_sizes_by_t0: dict,
    treated_sizes_by_t0: dict,
    pseudocount: float = 1.0,
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Perform differential enrichment analysis.

    For every T0 clone, computes fold-changes relative to T0 for both
    Control and Treated, then tests for differential enrichment using
    Fisher's exact test with FDR correction.

    Parameters
    ----------
    t0_sizes : dict
        {t0_clone_id: UMI_count}
    control_sizes_by_t0 : dict
        {t0_clone_id: control_UMI_count}
    treated_sizes_by_t0 : dict
        {t0_clone_id: treated_UMI_count}
    pseudocount : float
        Pseudocount added before fold-change calculation.
    fdr_threshold : float
        FDR q-value threshold for significance.
    log2fc_threshold : float
        Absolute Log2 fold-change ratio threshold.

    Returns
    -------
    pd.DataFrame
        Results table with columns defined in RESULT_COLUMNS.
    """
    all_t0_clones = sorted(t0_sizes.keys())
    n_clones = len(all_t0_clones)
    print(f"  Analysing {n_clones} T0 clones ...")

    records = []
    for cid in all_t0_clones:
        t0_sz = t0_sizes[cid]
        ctrl_sz = control_sizes_by_t0.get(cid, 0)
        treat_sz = treated_sizes_by_t0.get(cid, 0)

        fc_ctrl = (ctrl_sz + pseudocount) / (t0_sz + pseudocount)
        fc_treat = (treat_sz + pseudocount) / (t0_sz + pseudocount)

        # Avoid division by zero (both fc values are > 0 due to pseudocount)
        log2_fc_ratio = np.log2(fc_treat / fc_ctrl)

        # Fisher's exact test contingency table
        # [[T0→Control, T1_Control], [T0→Treated, T1_Treated]]
        # Use integer counts; add pseudocount as 1 to avoid zero cells
        table = [
            [int(t0_sz) + 1, int(ctrl_sz) + 1],
            [int(t0_sz) + 1, int(treat_sz) + 1],
        ]
        try:
            _, p_val = fisher_exact(table, alternative="two-sided")
        except Exception:
            p_val = 1.0

        records.append({
            "Clone_ID": cid,
            "T0_size": t0_sz,
            "Control_size": ctrl_sz,
            "Treated_size": treat_sz,
            "FC_control": fc_ctrl,
            "FC_treated": fc_treat,
            "Log2_FC_ratio": log2_fc_ratio,
            "p_value": p_val,
        })

    df = pd.DataFrame(records)

    # FDR correction (Benjamini-Hochberg)
    if len(df) > 0:
        reject, q_values, _, _ = multipletests(
            df["p_value"].values, alpha=fdr_threshold, method="fdr_bh"
        )
        df["q_value"] = q_values
    else:
        df["q_value"] = []

    # Classification
    def classify(row):
        if row["q_value"] < fdr_threshold and row["Log2_FC_ratio"] > log2fc_threshold:
            return "Enriched"
        elif row["q_value"] < fdr_threshold and row["Log2_FC_ratio"] < -log2fc_threshold:
            return "Depleted"
        else:
            return "Neutral"

    df["Status"] = df.apply(classify, axis=1)

    # Sort by q_value then absolute log2FC
    df = df.sort_values(["q_value", "p_value"], ascending=True).reset_index(drop=True)

    return df[RESULT_COLUMNS]


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_results(df: pd.DataFrame, output_dir: str):
    """Write all output files to *output_dir*."""
    os.makedirs(output_dir, exist_ok=True)

    # Main results
    df.to_csv(
        os.path.join(output_dir, "differential_enrichment_results.txt"),
        sep="\t", index=False,
    )

    # Subsets
    enriched = df[df["Status"] == "Enriched"]
    depleted = df[df["Status"] == "Depleted"]
    neutral = df[df["Status"] == "Neutral"]

    enriched.to_csv(os.path.join(output_dir, "enriched_clones.txt"), sep="\t", index=False)
    depleted.to_csv(os.path.join(output_dir, "depleted_clones.txt"), sep="\t", index=False)
    neutral.to_csv(os.path.join(output_dir, "neutral_clones.txt"), sep="\t", index=False)

    # Summary statistics
    n_total = len(df)
    n_enriched = len(enriched)
    n_depleted = len(depleted)
    n_neutral = len(neutral)

    lines = [
        "=" * 60,
        "Differential Enrichment Analysis — Summary",
        "=" * 60,
        f"Total T0 clones analysed:    {n_total}",
        f"Enriched clones:             {n_enriched} ({100 * n_enriched / max(n_total, 1):.1f}%)",
        f"Depleted clones:             {n_depleted} ({100 * n_depleted / max(n_total, 1):.1f}%)",
        f"Neutral clones:              {n_neutral} ({100 * n_neutral / max(n_total, 1):.1f}%)",
        "",
        "--- Fold-change statistics ---",
    ]

    if n_total > 0:
        fc = df["Log2_FC_ratio"]
        lines += [
            f"Log2_FC_ratio  mean:    {fc.mean():.4f}",
            f"Log2_FC_ratio  median:  {fc.median():.4f}",
            f"Log2_FC_ratio  std:     {fc.std():.4f}",
            f"Log2_FC_ratio  min:     {fc.min():.4f}",
            f"Log2_FC_ratio  max:     {fc.max():.4f}",
        ]

    lines += [
        "",
        "--- Dropout statistics ---",
        f"Clones absent in Control:  {(df['Control_size'] == 0).sum()}",
        f"Clones absent in Treated:  {(df['Treated_size'] == 0).sum()}",
        f"Clones absent in both:     {((df['Control_size'] == 0) & (df['Treated_size'] == 0)).sum()}",
    ]

    if n_enriched > 0:
        lines += [
            "",
            "--- Top 10 enriched clones ---",
        ]
        top_enriched = enriched.sort_values("Log2_FC_ratio", ascending=False).head(10)
        for _, r in top_enriched.iterrows():
            lines.append(
                f"  Clone {int(r['Clone_ID']):>6d}  Log2FC_ratio={r['Log2_FC_ratio']:+.3f}  "
                f"q={r['q_value']:.2e}  T0={int(r['T0_size'])}  Ctrl={int(r['Control_size'])}  "
                f"Treat={int(r['Treated_size'])}"
            )

    if n_depleted > 0:
        lines += [
            "",
            "--- Top 10 depleted clones ---",
        ]
        top_depleted = depleted.sort_values("Log2_FC_ratio", ascending=True).head(10)
        for _, r in top_depleted.iterrows():
            lines.append(
                f"  Clone {int(r['Clone_ID']):>6d}  Log2FC_ratio={r['Log2_FC_ratio']:+.3f}  "
                f"q={r['q_value']:.2e}  T0={int(r['T0_size'])}  Ctrl={int(r['Control_size'])}  "
                f"Treat={int(r['Treated_size'])}"
            )

    summary_path = os.path.join(output_dir, "summary_statistics.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"  Results written to {output_dir}")
    print(f"    Enriched: {n_enriched}  |  Depleted: {n_depleted}  |  Neutral: {n_neutral}")


# ---------------------------------------------------------------------------
# Visualisations
# ---------------------------------------------------------------------------

def _setup_plot_style():
    """Set consistent, publication-quality matplotlib style."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })


def _status_palette():
    return {"Enriched": "#d62728", "Depleted": "#1f77b4", "Neutral": "#999999"}


def generate_plots(
    df: pd.DataFrame,
    output_dir: str,
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
):
    """Generate all publication-quality plots."""
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    if len(df) == 0:
        print("  WARNING: No data to plot — skipping visualisations.")
        return

    _setup_plot_style()
    palette = _status_palette()

    # 1. Volcano plot
    _plot_volcano(df, plot_dir, palette, fdr_threshold, log2fc_threshold)

    # 2. MA plot
    _plot_ma(df, plot_dir, palette)

    # 3. Clone trajectory plot
    _plot_trajectory(df, plot_dir, palette)

    # 4. Enrichment distribution
    _plot_distribution(df, plot_dir)

    # 5. Top clones heatmap
    _plot_heatmap(df, plot_dir)

    # 6. Summary bar plot
    _plot_summary_bar(df, plot_dir, palette)

    print(f"  Plots saved to {plot_dir}")


def _plot_volcano(df, plot_dir, palette, fdr_thresh, fc_thresh):
    """Volcano plot: Log2_FC_ratio vs -log10(q_value)."""
    fig, ax = plt.subplots(figsize=(8, 6))

    neg_log_q = -np.log10(df["q_value"].clip(lower=1e-300))

    for status in ["Neutral", "Depleted", "Enriched"]:
        mask = df["Status"] == status
        ax.scatter(
            df.loc[mask, "Log2_FC_ratio"],
            neg_log_q[mask],
            c=palette[status],
            label=status,
            s=12, alpha=0.6, edgecolors="none",
        )

    # Threshold lines
    ax.axhline(-np.log10(fdr_thresh), ls="--", lw=0.8, color="grey", alpha=0.7)
    ax.axvline(fc_thresh, ls="--", lw=0.8, color="grey", alpha=0.7)
    ax.axvline(-fc_thresh, ls="--", lw=0.8, color="grey", alpha=0.7)

    ax.set_xlabel("Log₂ Fold-Change Ratio (Treated/Control)")
    ax.set_ylabel("-log₁₀(q-value)")
    ax.set_title("Volcano Plot — Differential Enrichment")
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "volcano_plot.png"))
    plt.close(fig)


def _plot_ma(df, plot_dir, palette):
    """MA plot: Log2(mean abundance) vs Log2_FC_ratio."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Mean abundance = average of Control and Treated sizes (+ 1 to avoid log(0))
    mean_abundance = np.log2(
        ((df["Control_size"] + df["Treated_size"]) / 2).clip(lower=0.5)
    )

    for status in ["Neutral", "Depleted", "Enriched"]:
        mask = df["Status"] == status
        ax.scatter(
            mean_abundance[mask],
            df.loc[mask, "Log2_FC_ratio"],
            c=palette[status],
            label=status,
            s=12, alpha=0.6, edgecolors="none",
        )

    ax.axhline(0, ls="-", lw=0.5, color="black", alpha=0.5)
    ax.set_xlabel("Log₂ Mean Abundance (Control + Treated)/2")
    ax.set_ylabel("Log₂ Fold-Change Ratio (Treated/Control)")
    ax.set_title("MA Plot — Enrichment vs Abundance")
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "ma_plot.png"))
    plt.close(fig)


def _plot_trajectory(df, plot_dir, palette):
    """Scatter of FC_control vs FC_treated (log2)."""
    fig, ax = plt.subplots(figsize=(7, 7))

    log2_fc_ctrl = np.log2(df["FC_control"].clip(lower=1e-10))
    log2_fc_treat = np.log2(df["FC_treated"].clip(lower=1e-10))

    for status in ["Neutral", "Depleted", "Enriched"]:
        mask = df["Status"] == status
        ax.scatter(
            log2_fc_ctrl[mask], log2_fc_treat[mask],
            c=palette[status], label=status,
            s=12, alpha=0.6, edgecolors="none",
        )

    # Diagonal
    lims = [
        min(log2_fc_ctrl.min(), log2_fc_treat.min()) - 0.5,
        max(log2_fc_ctrl.max(), log2_fc_treat.max()) + 0.5,
    ]
    ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel("Log₂ FC Control (T1_ctrl / T0)")
    ax.set_ylabel("Log₂ FC Treated (T1_treat / T0)")
    ax.set_title("Clone Trajectory — Control vs Treated Fold-Change")
    ax.legend(loc="upper left", framealpha=0.9)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "clone_trajectory_plot.png"))
    plt.close(fig)


def _plot_distribution(df, plot_dir):
    """Histogram of Log2_FC_ratio values."""
    fig, ax = plt.subplots(figsize=(8, 5))

    values = df["Log2_FC_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    ax.hist(values, bins=80, color="#5a9bd5", edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.axvline(0, ls="-", lw=1, color="black", alpha=0.5)
    ax.axvline(1.0, ls="--", lw=0.8, color="red", alpha=0.7, label="Enrichment threshold")
    ax.axvline(-1.0, ls="--", lw=0.8, color="blue", alpha=0.7, label="Depletion threshold")

    ax.set_xlabel("Log₂ Fold-Change Ratio (Treated/Control)")
    ax.set_ylabel("Number of Clones")
    ax.set_title("Distribution of Log₂ FC Ratios")
    ax.legend(framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "enrichment_distribution.png"))
    plt.close(fig)


def _plot_heatmap(df, plot_dir, n_top: int = 20):
    """Heatmap of top enriched and top depleted clones."""
    enriched = df[df["Status"] == "Enriched"].sort_values("Log2_FC_ratio", ascending=False).head(n_top)
    depleted = df[df["Status"] == "Depleted"].sort_values("Log2_FC_ratio", ascending=True).head(n_top)
    subset = pd.concat([enriched, depleted])

    if len(subset) == 0:
        return

    # Prepare heatmap data (Log2 sizes + 1)
    heat_data = subset[["Clone_ID", "T0_size", "Control_size", "Treated_size"]].copy()
    heat_data = heat_data.set_index("Clone_ID")
    heat_data.columns = ["T0", "Control", "Treated"]
    heat_data = np.log2(heat_data + 1)

    fig, ax = plt.subplots(figsize=(6, max(4, len(subset) * 0.3 + 1)))
    sns.heatmap(
        heat_data, ax=ax, cmap="YlOrRd", annot=False,
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Log₂(UMI count + 1)"},
    )
    ax.set_title(f"Top {n_top} Enriched & Depleted Clones")
    ax.set_ylabel("Clone ID")
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "top_clones_heatmap.png"))
    plt.close(fig)


def _plot_summary_bar(df, plot_dir, palette):
    """Bar chart of clone classification counts."""
    counts = df["Status"].value_counts()
    categories = ["Enriched", "Neutral", "Depleted"]
    values = [counts.get(c, 0) for c in categories]
    colors = [palette[c] for c in categories]

    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(categories, values, color=colors, edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.02,
            str(val), ha="center", va="bottom", fontweight="bold", fontsize=10,
        )
    ax.set_ylabel("Number of Clones")
    ax.set_title("Clone Classification Summary")
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "summary_bar_plot.png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Differential Enrichment Analysis for RewindPipeline\n"
            "===================================================\n"
            "Identifies clones significantly enriched or depleted by\n"
            "drug selection, controlling for random drift in the\n"
            "Control population.\n"
        ),
    )

    inp = parser.add_argument_group("Input files")
    inp.add_argument(
        "--t0-clones", required=True,
        help="T0 cloneID_cloneBarcode.tsv file (from Step 10)",
    )
    inp.add_argument(
        "--control-clones", required=True,
        help="T1 Control cloneID_cloneBarcode.tsv file",
    )
    inp.add_argument(
        "--treated-clones", required=True,
        help="T1 Treated cloneID_cloneBarcode.tsv file",
    )
    inp.add_argument(
        "--control-matches", required=True,
        help="Control best_matches.txt from match_timepoints.py",
    )
    inp.add_argument(
        "--treated-matches", required=True,
        help="Treated best_matches.txt from match_timepoints.py",
    )

    # Optional UMI count files for clone sizing
    inp.add_argument(
        "--t0-filtered", default=None,
        help="T0 rewind_filtered.tsv — used for UMI-based clone sizes. "
             "If not provided, clone sizes are estimated from the number "
             "of lineage barcodes per clone.",
    )
    inp.add_argument(
        "--control-filtered", default=None,
        help="Control rewind_filtered.tsv for UMI-based clone sizes.",
    )
    inp.add_argument(
        "--treated-filtered", default=None,
        help="Treated rewind_filtered.tsv for UMI-based clone sizes.",
    )

    out = parser.add_argument_group("Output")
    out.add_argument(
        "--output-dir", default="enrichment_results",
        help="Output directory (default: enrichment_results)",
    )

    stat = parser.add_argument_group("Statistical parameters")
    stat.add_argument(
        "--pseudocount", type=float, default=1.0,
        help="Pseudocount for fold-change calculation (default: 1)",
    )
    stat.add_argument(
        "--fdr-threshold", type=float, default=0.05,
        help="FDR q-value threshold for significance (default: 0.05)",
    )
    stat.add_argument(
        "--log2fc-threshold", type=float, default=1.0,
        help="Log2 fold-change ratio threshold (default: 1.0)",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print("=" * 60)
    print("  Differential Enrichment Analysis")
    print("=" * 60)

    # ---- 1. Load clone barcode data ----
    print("\n[1/5] Loading clone barcode data ...")
    t0_barcodes = parse_clone_barcode_file(args.t0_clones)
    ctrl_barcodes = parse_clone_barcode_file(args.control_clones)
    treat_barcodes = parse_clone_barcode_file(args.treated_clones)
    print(f"  T0 clones:      {len(t0_barcodes)}")
    print(f"  Control clones:  {len(ctrl_barcodes)}")
    print(f"  Treated clones:  {len(treat_barcodes)}")

    # ---- 2. Compute clone sizes ----
    print("\n[2/5] Computing clone sizes ...")

    if args.t0_filtered:
        t0_bc_umis = load_clone_sizes_from_filtered(args.t0_filtered)
        t0_sizes = compute_clone_sizes(t0_barcodes, t0_bc_umis)
        print(f"  T0 sizes: UMI-based (from {args.t0_filtered})")
    else:
        # Fallback: use number of barcodes as proxy for clone size
        t0_sizes = {cid: len(bcs) for cid, bcs in t0_barcodes.items()}
        print("  T0 sizes: barcode-count proxy (no --t0-filtered provided)")

    if args.control_filtered:
        ctrl_bc_umis = load_clone_sizes_from_filtered(args.control_filtered)
        ctrl_sizes = compute_clone_sizes(ctrl_barcodes, ctrl_bc_umis)
        print(f"  Control sizes: UMI-based (from {args.control_filtered})")
    else:
        ctrl_sizes = {cid: len(bcs) for cid, bcs in ctrl_barcodes.items()}
        print("  Control sizes: barcode-count proxy")

    if args.treated_filtered:
        treat_bc_umis = load_clone_sizes_from_filtered(args.treated_filtered)
        treat_sizes = compute_clone_sizes(treat_barcodes, treat_bc_umis)
        print(f"  Treated sizes: UMI-based (from {args.treated_filtered})")
    else:
        treat_sizes = {cid: len(bcs) for cid, bcs in treat_barcodes.items()}
        print("  Treated sizes: barcode-count proxy")

    # ---- 3. Load matching results and build T0-indexed maps ----
    print("\n[3/5] Loading timepoint matching results ...")
    ctrl_match = load_best_matches(args.control_matches)
    treat_match = load_best_matches(args.treated_matches)
    print(f"  Control matches: {len(ctrl_match)} T1→T0 links")
    print(f"  Treated matches: {len(treat_match)} T1→T0 links")

    # Aggregate T1 sizes by T0 clone ID
    ctrl_sizes_by_t0 = build_t0_to_t1_size_map(ctrl_sizes, ctrl_match)
    treat_sizes_by_t0 = build_t0_to_t1_size_map(treat_sizes, treat_match)
    print(f"  T0 clones with Control data: {len(ctrl_sizes_by_t0)}")
    print(f"  T0 clones with Treated data: {len(treat_sizes_by_t0)}")

    # ---- 4. Run enrichment analysis ----
    print("\n[4/5] Running differential enrichment analysis ...")
    results = run_enrichment_analysis(
        t0_sizes=t0_sizes,
        control_sizes_by_t0=ctrl_sizes_by_t0,
        treated_sizes_by_t0=treat_sizes_by_t0,
        pseudocount=args.pseudocount,
        fdr_threshold=args.fdr_threshold,
        log2fc_threshold=args.log2fc_threshold,
    )

    # ---- 5. Write results and generate plots ----
    print("\n[5/5] Writing results and generating plots ...")
    write_results(results, args.output_dir)
    generate_plots(
        results, args.output_dir,
        fdr_threshold=args.fdr_threshold,
        log2fc_threshold=args.log2fc_threshold,
    )

    print(f"\n{'='*60}")
    print("  Differential Enrichment Analysis complete!")
    print(f"  Results in: {os.path.abspath(args.output_dir)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
