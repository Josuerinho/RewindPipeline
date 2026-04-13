#!/usr/bin/env python3
"""
Timepoint Matching for RewindPipeline
======================================
Matches clones between two timepoints (T0 → T1) using fuzzy barcode
matching with configurable Hamming distance and Fuzzy Jaccard similarity.

Algorithm overview
------------------
1. Load clone → barcode mappings from both timepoints (the
   ``*_cloneID_cloneBarcode.tsv`` output of Step 10).
2. Build a lookup structure mapping individual lineage barcodes → clone IDs.
3. For each T1 clone, find candidate T0 clones via barcode overlap:
   a. For each T1 barcode, find all T0 barcodes within the Hamming distance
      threshold (fuzzy matching).
   b. Collect all T0 clones that contain any of those matched T0 barcodes.
   c. Compute a Fuzzy Jaccard similarity for each T0 candidate:
      score = |matched_barcodes| / |T1_barcodes ∪ T0_barcodes|
      Optionally weighted by UMI counts.
4. Rank candidate matches and apply a similarity threshold.

Outputs
-------
- best_matches.txt       – one best T0 match per T1 clone
- ranked_matches.txt     – top-N T0 matches per T1 clone
- unmatched_clones.txt   – T1 clones below similarity threshold
- matching_statistics.txt – summary statistics
- plots/                 – QC visualisation plots

Usage
-----
    python match_timepoints.py \
        --t0-clone-barcode T0_cloneID_cloneBarcode.tsv \
        --t1-clone-barcode T1_cloneID_cloneBarcode.tsv \
        --output-dir timepoint_matching_results

    # Optionally supply filtered TSVs for UMI-weighted scoring:
    python match_timepoints.py \
        --t0-clone-barcode T0_cloneID_cloneBarcode.tsv \
        --t1-clone-barcode T1_cloneID_cloneBarcode.tsv \
        --t0-filtered T0_rewind_filtered.tsv \
        --t1-filtered T1_rewind_filtered.tsv \
        --weight-by-umi \
        --output-dir timepoint_matching_results
"""

import argparse
import os
import sys
import time
import warnings
from collections import defaultdict
from itertools import product
from multiprocessing import Pool, cpu_count

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def hamming_distance(s1: str, s2: str) -> int:
    """Return the Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))  # treat unequal lengths as max dist
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def parse_clone_barcode_file(path: str) -> dict:
    """
    Parse a cloneID_cloneBarcode.tsv file.

    Returns
    -------
    dict
        {clone_id (int): set of individual lineage barcode strings}
    """
    df = pd.read_csv(path, sep="\t")
    clone_to_barcodes: dict[int, set] = {}
    for _, row in df.iterrows():
        cid = int(row["cloneID"])
        barcodes = set(str(row["clone_barcode"]).split("_"))
        clone_to_barcodes[cid] = barcodes
    return clone_to_barcodes


def load_umi_weights(filtered_path: str) -> dict:
    """
    Load UMI counts from a rewind_filtered.tsv file.

    Returns
    -------
    dict
        {lineage_barcode: total_umi_count}  (summed across all cells)
    """
    df = pd.read_csv(filtered_path, sep="\t")
    umi_col = "umi_count" if "umi_count" in df.columns else df.columns[2]
    bc_col = "corrected_lineage_barcode" if "corrected_lineage_barcode" in df.columns else df.columns[1]
    weights = df.groupby(bc_col)[umi_col].sum().to_dict()
    return weights


def build_barcode_index(clone_to_barcodes: dict) -> dict:
    """
    Build barcode → set of clone IDs index.

    Returns
    -------
    dict
        {barcode_str: set of clone_ids}
    """
    idx: dict[str, set] = defaultdict(set)
    for cid, bcs in clone_to_barcodes.items():
        for bc in bcs:
            idx[bc].add(cid)
    return dict(idx)


def find_fuzzy_matches(query_bc: str, ref_barcodes: list, threshold: int) -> list:
    """
    Find all reference barcodes within *threshold* Hamming distance of *query_bc*.

    Returns list of (ref_barcode, hamming_dist) tuples.
    """
    matches = []
    for ref_bc in ref_barcodes:
        d = hamming_distance(query_bc, ref_bc)
        if d <= threshold:
            matches.append((ref_bc, d))
    return matches


# ---------------------------------------------------------------------------
# Core matching logic
# ---------------------------------------------------------------------------

def _match_single_t1_clone(args):
    """Worker function for parallel matching of a single T1 clone."""
    (
        t1_cid,
        t1_barcodes,
        t0_barcode_list,
        t0_barcode_to_clones,
        t0_clone_to_barcodes,
        hamming_threshold,
        top_n,
        weight_by_umi,
        t0_umi_weights,
        t1_umi_weights,
    ) = args

    # Step 1: For each T1 barcode, find fuzzy-matching T0 barcodes
    t1_to_t0_matches: dict[str, list] = {}  # t1_bc -> [(t0_bc, dist), ...]
    candidate_t0_clones: set = set()

    for t1_bc in t1_barcodes:
        matches = find_fuzzy_matches(t1_bc, t0_barcode_list, hamming_threshold)
        t1_to_t0_matches[t1_bc] = matches
        for t0_bc, _ in matches:
            candidate_t0_clones.update(t0_barcode_to_clones.get(t0_bc, set()))

    # Step 2: Score each candidate T0 clone
    results = []
    for t0_cid in candidate_t0_clones:
        t0_barcodes = t0_clone_to_barcodes[t0_cid]

        # Find matched barcode pairs
        matched_t1_bcs = set()
        matched_t0_bcs = set()
        for t1_bc in t1_barcodes:
            for t0_bc, dist in t1_to_t0_matches.get(t1_bc, []):
                if t0_bc in t0_barcodes:
                    matched_t1_bcs.add(t1_bc)
                    matched_t0_bcs.add(t0_bc)

        num_matched = len(matched_t1_bcs)
        if num_matched == 0:
            continue

        # Fuzzy Jaccard similarity
        union_size = len(t1_barcodes) + len(t0_barcodes) - num_matched
        if union_size == 0:
            continue
        similarity = num_matched / union_size

        # Optional UMI weighting
        umi_weighted_score = None
        if weight_by_umi and t0_umi_weights and t1_umi_weights:
            matched_umi_sum = sum(
                t1_umi_weights.get(bc, 1) for bc in matched_t1_bcs
            ) + sum(
                t0_umi_weights.get(bc, 1) for bc in matched_t0_bcs
            )
            total_umi_sum = sum(
                t1_umi_weights.get(bc, 1) for bc in t1_barcodes
            ) + sum(
                t0_umi_weights.get(bc, 1) for bc in t0_barcodes
            )
            umi_weighted_score = matched_umi_sum / total_umi_sum if total_umi_sum > 0 else 0.0

        results.append({
            "T1_clone_id": t1_cid,
            "T0_clone_id": t0_cid,
            "similarity_score": similarity,
            "umi_weighted_score": umi_weighted_score,
            "num_matched_barcodes": num_matched,
            "num_T1_barcodes": len(t1_barcodes),
            "num_T0_barcodes": len(t0_barcodes),
        })

    # Sort by similarity descending, keep top N
    results.sort(key=lambda x: x["similarity_score"], reverse=True)
    top_results = results[:top_n]

    # Always return t1_cid so caller can track zero-match clones
    return (t1_cid, top_results)


def match_timepoints(
    t0_clone_to_barcodes: dict,
    t1_clone_to_barcodes: dict,
    hamming_threshold: int = 2,
    similarity_threshold: float = 0.8,
    top_n: int = 5,
    weight_by_umi: bool = False,
    t0_umi_weights: dict | None = None,
    t1_umi_weights: dict | None = None,
    n_jobs: int = 1,
) -> tuple:
    """
    Match T1 clones to T0 clones using fuzzy barcode matching.

    Returns
    -------
    ranked_matches : list of list of dict
        For each T1 clone, a list of top-N candidate matches.
    """
    # Build index structures
    t0_barcode_to_clones = build_barcode_index(t0_clone_to_barcodes)
    t0_barcode_list = list(t0_barcode_to_clones.keys())

    # Prepare worker arguments
    worker_args = [
        (
            t1_cid,
            t1_barcodes,
            t0_barcode_list,
            t0_barcode_to_clones,
            t0_clone_to_barcodes,
            hamming_threshold,
            top_n,
            weight_by_umi,
            t0_umi_weights,
            t1_umi_weights,
        )
        for t1_cid, t1_barcodes in t1_clone_to_barcodes.items()
    ]

    total = len(worker_args)
    print(f"  Matching {total} T1 clones against {len(t0_clone_to_barcodes)} T0 clones ...")
    print(f"  T0 unique barcodes: {len(t0_barcode_list)}")
    print(f"  Hamming threshold: {hamming_threshold}")
    print(f"  Similarity threshold: {similarity_threshold}")
    print(f"  Using {n_jobs} worker(s)")

    t_start = time.time()

    # all_ranked: list of (t1_cid, [match_dicts])
    if n_jobs > 1:
        with Pool(n_jobs) as pool:
            all_ranked = []
            for i, result in enumerate(pool.imap(_match_single_t1_clone, worker_args, chunksize=50)):
                all_ranked.append(result)
                if (i + 1) % 500 == 0 or (i + 1) == total:
                    elapsed = time.time() - t_start
                    pct = 100 * (i + 1) / total
                    print(f"    Progress: {i+1}/{total} ({pct:.1f}%) — {elapsed:.1f}s", flush=True)
    else:
        all_ranked = []
        for i, wa in enumerate(worker_args):
            result = _match_single_t1_clone(wa)
            all_ranked.append(result)
            if (i + 1) % 500 == 0 or (i + 1) == total:
                elapsed = time.time() - t_start
                pct = 100 * (i + 1) / total
                print(f"    Progress: {i+1}/{total} ({pct:.1f}%) — {elapsed:.1f}s", flush=True)

    elapsed = time.time() - t_start
    print(f"  Matching complete in {elapsed:.1f}s")

    return all_ranked


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_results(
    all_ranked: list,
    similarity_threshold: float,
    output_dir: str,
    weight_by_umi: bool = False,
) -> tuple:
    """
    Write result files and return (best_df, ranked_df, unmatched_df) for
    downstream plotting.
    """
    os.makedirs(output_dir, exist_ok=True)

    # all_ranked is list of (t1_cid, [match_dicts])
    # Flatten all ranked matches
    flat = []
    zero_candidate_ids = []
    for t1_cid, matches in all_ranked:
        if matches:
            flat.extend(matches)
        else:
            zero_candidate_ids.append(t1_cid)

    if not flat:
        print("  WARNING: No matches found at all!")
        empty_df = pd.DataFrame()
        return empty_df, empty_df, empty_df

    ranked_df = pd.DataFrame(flat)

    # --- ranked_matches.txt ---
    cols = ["T1_clone_id", "T0_clone_id", "similarity_score",
            "num_matched_barcodes", "num_T1_barcodes", "num_T0_barcodes"]
    if weight_by_umi:
        cols.insert(3, "umi_weighted_score")
    ranked_df[cols].to_csv(
        os.path.join(output_dir, "ranked_matches.txt"), sep="\t", index=False
    )

    # --- best_matches.txt ---
    best_idx = ranked_df.groupby("T1_clone_id")["similarity_score"].idxmax()
    best_df = ranked_df.loc[best_idx].copy()
    best_df[cols].to_csv(
        os.path.join(output_dir, "best_matches.txt"), sep="\t", index=False
    )

    # --- unmatched_clones.txt ---
    # Clones with candidates but below threshold
    below_threshold = best_df[best_df["similarity_score"] < similarity_threshold][
        ["T1_clone_id", "T0_clone_id", "similarity_score"]
    ].copy()
    # Also include zero-candidate clones
    zero_rows = pd.DataFrame({
        "T1_clone_id": zero_candidate_ids,
        "T0_clone_id": [pd.NA] * len(zero_candidate_ids),
        "similarity_score": [0.0] * len(zero_candidate_ids),
    })
    unmatched_df = pd.concat([below_threshold, zero_rows], ignore_index=True)
    unmatched_df.to_csv(
        os.path.join(output_dir, "unmatched_clones.txt"), sep="\t", index=False
    )

    # --- matching_statistics.txt ---
    n_t1_total = len(all_ranked)
    n_t1_with_candidates = sum(1 for _, m in all_ranked if m)
    matched_above = set(best_df[best_df["similarity_score"] >= similarity_threshold]["T1_clone_id"])
    n_t1_matched = len(matched_above)
    n_t1_unmatched = n_t1_total - n_t1_matched

    scores = best_df["similarity_score"].values
    stats_lines = [
        f"Total T1 clones:               {n_t1_total}",
        f"T1 clones with candidates:     {n_t1_with_candidates}",
        f"T1 clones matched (>={similarity_threshold}): {n_t1_matched}",
        f"T1 clones unmatched:           {n_t1_unmatched}",
        f"Match rate:                    {100 * n_t1_matched / max(n_t1_total, 1):.1f}%",
        f"Mean similarity (best match):  {np.mean(scores):.4f}",
        f"Median similarity:             {np.median(scores):.4f}",
        f"Std similarity:                {np.std(scores):.4f}",
        f"Min similarity:                {np.min(scores):.4f}",
        f"Max similarity:                {np.max(scores):.4f}",
    ]
    if weight_by_umi:
        umi_scores = best_df["umi_weighted_score"].dropna().values
        if len(umi_scores):
            stats_lines.append(f"Mean UMI-weighted similarity:  {np.mean(umi_scores):.4f}")
            stats_lines.append(f"Median UMI-weighted similarity:{np.median(umi_scores):.4f}")

    stats_text = "\n".join(stats_lines)
    with open(os.path.join(output_dir, "matching_statistics.txt"), "w") as fh:
        fh.write(stats_text + "\n")
    print(f"\n  === Matching Statistics ===\n{stats_text}\n")

    return best_df, ranked_df, below_threshold


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def generate_plots(
    best_df: pd.DataFrame,
    ranked_df: pd.DataFrame,
    unmatched_df: pd.DataFrame,
    similarity_threshold: float,
    n_total_t1: int,
    output_dir: str,
):
    """Generate QC plots and save to *output_dir*/plots/."""
    if best_df.empty:
        print("  Skipping plots — no match data available.")
        return

    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    sns.set_style("whitegrid")

    # 1. Histogram of best-match similarity scores
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(best_df["similarity_score"], bins=50, edgecolor="black", alpha=0.75,
            color="steelblue")
    ax.axvline(similarity_threshold, color="red", linestyle="--", linewidth=1.5,
               label=f"Threshold = {similarity_threshold}")
    ax.set_xlabel("Best-match Fuzzy Jaccard Similarity")
    ax.set_ylabel("Number of T1 Clones")
    ax.set_title("Distribution of Best-Match Similarity Scores")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "similarity_histogram.png"), dpi=150)
    plt.close(fig)

    # 2. Scatter: T0 clone size vs T1 clone size (matched pairs)
    matched = best_df[best_df["similarity_score"] >= similarity_threshold].copy()
    if not matched.empty:
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(
            matched["num_T0_barcodes"],
            matched["num_T1_barcodes"],
            alpha=0.4,
            s=20,
            c=matched["similarity_score"],
            cmap="viridis",
            edgecolors="none",
        )
        cbar = fig.colorbar(ax.collections[0], ax=ax, label="Similarity score")
        ax.set_xlabel("T0 Clone Size (# barcodes)")
        ax.set_ylabel("T1 Clone Size (# barcodes)")
        ax.set_title("Matched Clone Sizes: T0 vs T1")
        # Add y=x reference line
        lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, lim], [0, lim], "k--", alpha=0.3, linewidth=1)
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, "clone_size_scatter.png"), dpi=150)
        plt.close(fig)

    # 3. Bar plot: matched vs unmatched
    n_matched = len(matched)
    n_unmatched = n_total_t1 - n_matched
    fig, ax = plt.subplots(figsize=(5, 5))
    bars = ax.bar(
        ["Matched\n(≥ threshold)", "Unmatched\n(< threshold)"],
        [n_matched, n_unmatched],
        color=["#2ca02c", "#d62728"],
        edgecolor="black",
    )
    for bar, val in zip(bars, [n_matched, n_unmatched]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(n_matched, n_unmatched) * 0.02,
            str(val),
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )
    ax.set_ylabel("Number of T1 Clones")
    ax.set_title(f"Match Status (threshold = {similarity_threshold})")
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "matched_vs_unmatched.png"), dpi=150)
    plt.close(fig)

    # 4. Heatmap of top clone pairs (if manageable)
    TOP_HEATMAP = 30  # show top 30 pairs
    if len(matched) > 0:
        top_pairs = matched.nlargest(min(TOP_HEATMAP, len(matched)), "similarity_score")
        pivot_data = top_pairs.pivot_table(
            index="T1_clone_id",
            columns="T0_clone_id",
            values="similarity_score",
            aggfunc="first",
        ).fillna(0)
        if pivot_data.shape[0] > 1 and pivot_data.shape[1] > 1:
            fig, ax = plt.subplots(
                figsize=(max(8, pivot_data.shape[1] * 0.5), max(6, pivot_data.shape[0] * 0.4))
            )
            sns.heatmap(
                pivot_data,
                annot=True,
                fmt=".2f",
                cmap="YlOrRd",
                linewidths=0.5,
                ax=ax,
                vmin=0,
                vmax=1,
            )
            ax.set_title(f"Top {min(TOP_HEATMAP, len(top_pairs))} Matched Clone Pairs")
            ax.set_xlabel("T0 Clone ID")
            ax.set_ylabel("T1 Clone ID")
            fig.tight_layout()
            fig.savefig(os.path.join(plot_dir, "top_matches_heatmap.png"), dpi=150)
            plt.close(fig)

    # 5. Similarity score vs number of matched barcodes
    if not matched.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(
            matched["num_matched_barcodes"],
            matched["similarity_score"],
            alpha=0.4,
            s=20,
            color="steelblue",
        )
        ax.set_xlabel("Number of Matched Barcodes")
        ax.set_ylabel("Fuzzy Jaccard Similarity")
        ax.set_title("Similarity vs Matched Barcode Count")
        fig.tight_layout()
        fig.savefig(os.path.join(plot_dir, "similarity_vs_matched_count.png"), dpi=150)
        plt.close(fig)

    print(f"  Plots saved to {plot_dir}/")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Match clones between two timepoints using fuzzy barcode matching.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    io = parser.add_argument_group("Input / Output")
    io.add_argument(
        "--t0-clone-barcode", required=True,
        help="Path to T0 cloneID_cloneBarcode.tsv (Step 10 output)",
    )
    io.add_argument(
        "--t1-clone-barcode", required=True,
        help="Path to T1 cloneID_cloneBarcode.tsv (Step 10 output)",
    )
    io.add_argument(
        "--t0-filtered", default=None,
        help="(Optional) T0 rewind_filtered.tsv — needed for --weight-by-umi",
    )
    io.add_argument(
        "--t1-filtered", default=None,
        help="(Optional) T1 rewind_filtered.tsv — needed for --weight-by-umi",
    )
    io.add_argument(
        "--output-dir", default="timepoint_matching_results",
        help="Output directory (default: timepoint_matching_results)",
    )

    params = parser.add_argument_group("Matching parameters")
    params.add_argument(
        "--hamming-threshold", type=int, default=2,
        help="Maximum Hamming distance for fuzzy barcode matching (default: 2)",
    )
    params.add_argument(
        "--similarity-threshold", type=float, default=0.8,
        help="Minimum Fuzzy Jaccard similarity to consider a match (default: 0.8)",
    )
    params.add_argument(
        "--top-n-matches", type=int, default=5,
        help="Number of top-ranked T0 matches to keep per T1 clone (default: 5)",
    )
    params.add_argument(
        "--weight-by-umi", action="store_true",
        help="Weight similarity by UMI counts (requires --t0-filtered and --t1-filtered)",
    )

    perf = parser.add_argument_group("Performance")
    perf.add_argument(
        "--n-jobs", type=int, default=1,
        help="Number of parallel workers (default: 1; set to 0 for auto-detect)",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("  Timepoint Matching — RewindPipeline")
    print("=" * 60)

    # --- Load data ---
    print("\n[1/4] Loading clone barcode files ...")
    t0_clone_to_barcodes = parse_clone_barcode_file(args.t0_clone_barcode)
    t1_clone_to_barcodes = parse_clone_barcode_file(args.t1_clone_barcode)
    print(f"  T0: {len(t0_clone_to_barcodes)} clones")
    print(f"  T1: {len(t1_clone_to_barcodes)} clones")

    # UMI weights
    t0_umi_weights = None
    t1_umi_weights = None
    if args.weight_by_umi:
        if not args.t0_filtered or not args.t1_filtered:
            sys.exit("ERROR: --weight-by-umi requires --t0-filtered and --t1-filtered")
        print("  Loading UMI weights ...")
        t0_umi_weights = load_umi_weights(args.t0_filtered)
        t1_umi_weights = load_umi_weights(args.t1_filtered)

    # --- Matching ---
    print("\n[2/4] Running fuzzy clone matching ...")
    n_jobs = args.n_jobs if args.n_jobs > 0 else max(1, cpu_count() - 1)
    all_ranked = match_timepoints(
        t0_clone_to_barcodes=t0_clone_to_barcodes,
        t1_clone_to_barcodes=t1_clone_to_barcodes,
        hamming_threshold=args.hamming_threshold,
        similarity_threshold=args.similarity_threshold,
        top_n=args.top_n_matches,
        weight_by_umi=args.weight_by_umi,
        t0_umi_weights=t0_umi_weights,
        t1_umi_weights=t1_umi_weights,
        n_jobs=n_jobs,
    )

    # --- Write results ---
    print("\n[3/4] Writing output files ...")
    best_df, ranked_df, unmatched_df = write_results(
        all_ranked=all_ranked,
        similarity_threshold=args.similarity_threshold,
        output_dir=args.output_dir,
        weight_by_umi=args.weight_by_umi,
    )

    # --- Plots ---
    print("\n[4/4] Generating QC plots ...")
    generate_plots(
        best_df=best_df,
        ranked_df=ranked_df,
        unmatched_df=unmatched_df,
        similarity_threshold=args.similarity_threshold,
        n_total_t1=len(t1_clone_to_barcodes),
        output_dir=args.output_dir,
    )

    print("\nDone! Results in:", os.path.abspath(args.output_dir))


if __name__ == "__main__":
    main()
