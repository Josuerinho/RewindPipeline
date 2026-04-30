import argparse
import itertools
import multiprocessing as mp
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import fsolve
from scipy.stats import kstest, poisson


# Globals for multiprocessing workers
_WITHIN_CELL_DATA = None
_CROSS_CELL_DATA_1 = None
_CROSS_CELL_DATA_2 = None
_WITHIN_PARAMS = None
_CROSS_PARAMS = None

DETAILED_PAIRWISE_COLUMNS = [
    "cell_barcode_1",
    "cell_barcode_2",
    "source_condition_1",
    "source_condition_2",
    "is_clone",
    "num_matches",
    "num_exact_matches",
    "num_hamming_matches",
    "num_barcode_comparisons",
    "min_observed_distance",
    "best_hamming_distance",
    "mean_match_confidence",
    "matched_barcode_pairs",
]


def parse_args():
    """
    Parses command-line arguments for the rewind_barcode_filtering_and_clone_assignment script.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Process lineage barcode data to filter noise and assign clones using "
            "Hamming-distance matching. Supports both within-condition and "
            "cross-condition comparisons."
        )
    )
    parser.add_argument(
        "--input-file",
        type=str,
        required=True,
        help=(
            "Path to input TSV containing starcode-corrected lineage barcode "
            "address data."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help=(
            "Path four output files. If not provided, outputs will be written "
            "to the same directory as the input file with appropriate suffixes."
        ),
    )

    parser.add_argument(
        "--input-file-2",
        type=str,
        default=None,
        help=(
            "Optional second input TSV used for cross-condition clone assignment. "
            "Required when --cross-condition is enabled."
        ),
    )
    parser.add_argument(
        "--cross-condition",
        action="store_true",
        help="Enable cross-condition clone assignment between --input-file and --input-file-2.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=1,
        help=(
            "Number of top log-normalized counts for expected count calculation "
            "(default: 1)."
        ),
    )
    parser.add_argument(
        "--z-threshold",
        type=float,
        default=2,
        help="Z-score threshold for filtering out noise (default: 2).",
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=5,
        help="Minimum group size for robust Z-score estimation (default: 5).",
    )
    parser.add_argument(
        "--max-lineage-count",
        type=int,
        default=None,
        help=(
            "Maximum lineage count to include for clone assignment. Cells with "
            "lineage counts >= this value are excluded from clone assignment."
        ),
    )

    parser.add_argument(
        "--hamming-threshold",
        type=int,
        default=2,
        help="Maximum Hamming distance for a barcode match (default: 2).",
    )
    parser.add_argument(
        "--quick-mode",
        dest="quick_mode",
        action="store_true",
        default=True,
        help=(
            "Enable set-intersection optimization: if exact barcode overlap exists "
            "for a cell pair, skip full pairwise Hamming comparisons (default: enabled)."
        ),
    )
    parser.add_argument(
        "--no-quick-mode",
        dest="quick_mode",
        action="store_false",
        help="Disable set-intersection optimization and always scan barcode pairs.",
    )
    parser.add_argument(
        "--force-hamming-all",
        action="store_true",
        help=(
            "Force full pairwise Hamming comparisons even when exact overlaps are "
            "found in quick mode."
        ),
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help=(
            "Number of parallel worker processes for pairwise comparisons "
            "(default: 1; use 0 for auto-detect)."
        ),
    )

    return parser.parse_args()


def filter_noise(df, top_n, z_threshold, min_group_size, target_sum=1e4):
    """
    Filter out noisy lineage barcode molecule counts using normalization and Z-scores.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with at least cell_barcode, umi_count, corrected_lineage_barcode.
    top_n : int
        Number of top log-normalized counts used for expected-count estimation.
    z_threshold : float
        Absolute Z-score threshold for filtering.
    min_group_size : int
        Minimum subpopulation size for robust Z-score computation.
    target_sum : float, optional
        Target sum for per-cell UMI normalization, by default 1e4.

    Returns
    -------
    tuple
        (df_with_stats, filtered_df)
    """
    df = df.copy()

    # Normalize UMI counts for each cell to sum to target_sum
    df["normalized_umi_count"] = df.groupby("cell_barcode")["umi_count"].transform(
        lambda x: (x / x.sum()) * target_sum
    )

    # Log-transform the normalized counts
    df["log_normalized_umi_count"] = np.log1p(df["normalized_umi_count"])

    # Expected count = mean of top-N highest log-normalized counts per cell
    df["rank"] = df.groupby("cell_barcode")["log_normalized_umi_count"].rank(
        method="first", ascending=False
    )
    expected_counts = (
        df[df["rank"] <= top_n]
        .groupby("cell_barcode")["log_normalized_umi_count"]
        .mean()
        .reset_index()
        .rename(columns={"log_normalized_umi_count": "expected_count"})
    )
    df = pd.merge(df, expected_counts, on="cell_barcode")

    # Subpopulation-level std dev by number of unique lineages per cell
    lineage_count_per_cell = (
        df.groupby("cell_barcode")["corrected_lineage_barcode"]
        .nunique()
        .reset_index(name="lineage_count")
    )
    df = pd.merge(df, lineage_count_per_cell, on="cell_barcode")
    std_devs = (
        df.groupby("lineage_count")["log_normalized_umi_count"]
        .std()
        .reset_index(name="std_dev")
    )
    df = pd.merge(df, std_devs, on="lineage_count", how="left")

    # Identify small groups
    group_sizes = df.groupby("cell_barcode")["lineage_count"].transform("size")
    df["is_small_group"] = group_sizes < min_group_size

    def calculate_z_score(row):
        if row["is_small_group"]:
            return np.nan
        if row["std_dev"] == 0 and row["lineage_count"] == 1:
            return 0
        return (row["log_normalized_umi_count"] - row["expected_count"]) / row["std_dev"]

    df["z_score"] = df.apply(calculate_z_score, axis=1)

    filtered_df = df[~df["z_score"].isna() & (df["z_score"].abs() < z_threshold)].copy()
    filtered_df["lineage_count"] = filtered_df.groupby("cell_barcode")[
        "corrected_lineage_barcode"
    ].transform("nunique")

    return df, filtered_df


def test_truncated_poisson_distribution(total_df, column_name):
    """
    Test fit of lineage counts to a zero-truncated Poisson distribution.

    Parameters
    ----------
    total_df : pd.DataFrame
        Dataframe containing the input count column.
    column_name : str
        Name of the count column.

    Returns
    -------
    tuple
        (empirical_counts, expected_counts, empirical_cdf, expected_cdf, results_dict)
    """
    empirical_counts = total_df[column_name].value_counts().sort_index()
    empirical_cdf = np.cumsum(empirical_counts / empirical_counts.sum())

    observed_mean = total_df[column_name].mean()

    def solve_lambda(lambda_estimate):
        return lambda_estimate / (1 - np.exp(-lambda_estimate)) - observed_mean

    lambda_estimate = fsolve(solve_lambda, observed_mean)[0]

    truncation_correction = 1 - poisson.pmf(0, lambda_estimate)
    max_lineages = empirical_counts.index.max()
    poisson_dist = (
        poisson.pmf(k=np.arange(1, max_lineages + 1), mu=lambda_estimate)
        / truncation_correction
    )
    expected_counts = poisson_dist * empirical_counts.sum()
    expected_cdf = np.cumsum(expected_counts) / empirical_counts.sum()

    min_length = min(len(expected_cdf), len(empirical_cdf))
    expected_counts = expected_counts[:min_length]
    expected_cdf = expected_cdf[:min_length]
    empirical_counts = empirical_counts[:min_length]
    empirical_cdf = empirical_cdf[:min_length]

    ks_stat, p_value = kstest(empirical_cdf, expected_cdf)
    percent_with_lentivirus = (1 - np.exp(-lambda_estimate)) * 100

    results = {
        "Observed Mean": observed_mean,
        "Estimated Lambda": lambda_estimate,
        "K-S Statistic": ks_stat,
        "P-Value": p_value,
        "Percent with Lentivirus": percent_with_lentivirus,
    }
    return empirical_counts, expected_counts, empirical_cdf, expected_cdf, results


def calculate_hamming_distance(seq1, seq2):
    """
    Calculate Hamming distance between two sequences.

    For unequal-length inputs, uses max length as a conservative distance.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.

    Returns
    -------
    int
        Hamming distance.
    """
    if len(seq1) != len(seq2):
        return max(len(seq1), len(seq2))
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def calculate_match_confidence(umi1, umi2, rank1, rank2):
    """
    Calculate confidence score for a barcode match.

    Score combines UMI balance and barcode rank agreement within each cell.

    Parameters
    ----------
    umi1 : float
        UMI count of barcode in first cell.
    umi2 : float
        UMI count of barcode in second cell.
    rank1 : float
        Rank of barcode in first cell (1 = strongest barcode).
    rank2 : float
        Rank of barcode in second cell.

    Returns
    -------
    float
        Confidence score in [0, 1].
    """
    umi1 = float(max(0, umi1))
    umi2 = float(max(0, umi2))
    max_umi = max(umi1, umi2)
    umi_balance = (min(umi1, umi2) / max_umi) if max_umi > 0 else 0.0

    rank_delta = abs(float(rank1) - float(rank2))
    rank_agreement = 1.0 / (1.0 + rank_delta)

    confidence = 0.7 * umi_balance + 0.3 * rank_agreement
    return max(0.0, min(1.0, confidence))


def _prepare_cell_barcode_data(input_df, max_lineage_count=None):
    """Prepare per-cell barcode records and exact-match lookup sets."""
    if input_df.empty:
        return {}, {}

    work_df = input_df.copy()

    if max_lineage_count is not None:
        work_df = work_df[work_df["lineage_count"] < max_lineage_count].copy()

    if work_df.empty:
        return {}, {}

    # Rank barcodes within each cell by UMI abundance.
    work_df["barcode_rank"] = work_df.groupby("cell_barcode")["umi_count"].rank(
        method="first", ascending=False
    )

    grouped = {}
    for cell, cell_df in work_df.groupby("cell_barcode", sort=True):
        entries = [
            {
                "barcode": str(r.corrected_lineage_barcode),
                "umi_count": float(r.umi_count),
                "rank": float(r.barcode_rank),
            }
            for r in cell_df.itertuples(index=False)
        ]
        grouped[cell] = entries

    barcode_sets = {
        cell: {entry["barcode"] for entry in entries}
        for cell, entries in grouped.items()
    }

    return grouped, barcode_sets


def _serialize_matches(matches):
    """Serialize match dictionaries into a compact string for output files."""
    if not matches:
        return ""
    return ";".join(
        [
            (
                f"{m['barcode_1']}|{m['barcode_2']}|d:{m['hamming_distance']}"
                f"|c:{m['confidence']:.4f}|u1:{m['umi_1']:.0f}|u2:{m['umi_2']:.0f}"
                f"|r1:{m['rank_1']:.0f}|r2:{m['rank_2']:.0f}"
            )
            for m in matches
        ]
    )


def _compare_two_cells(cell1, cell2, entries1, entries2, set1, set2, hamming_threshold, quick_mode, force_hamming_all, source_1, source_2):
    """Compare lineage barcodes between two cells and return detailed match diagnostics."""
    if cell1 == cell2 and source_1 == source_2:
        return {
            "cell_barcode_1": cell1,
            "cell_barcode_2": cell2,
            "source_condition_1": source_1,
            "source_condition_2": source_2,
            "is_clone": False,
            "num_matches": 0,
            "num_exact_matches": 0,
            "num_hamming_matches": 0,
            "num_barcode_comparisons": 0,
            "min_observed_distance": np.nan,
            "best_hamming_distance": np.nan,
            "mean_match_confidence": np.nan,
            "matched_barcode_pairs": "",
        }

    matches = []
    num_barcode_comparisons = 0
    min_observed_distance = np.nan

    # Quick optimization: exact intersections can be accepted immediately.
    exact_intersection = set1.intersection(set2) if quick_mode else set()
    if exact_intersection:
        entry1_map = {e["barcode"]: e for e in entries1}
        entry2_map = {e["barcode"]: e for e in entries2}

        for bc in sorted(exact_intersection):
            e1 = entry1_map[bc]
            e2 = entry2_map[bc]
            conf = calculate_match_confidence(e1["umi_count"], e2["umi_count"], e1["rank"], e2["rank"])
            matches.append(
                {
                    "barcode_1": bc,
                    "barcode_2": bc,
                    "hamming_distance": 0,
                    "confidence": conf,
                    "umi_1": e1["umi_count"],
                    "umi_2": e2["umi_count"],
                    "rank_1": e1["rank"],
                    "rank_2": e2["rank"],
                }
            )
            num_barcode_comparisons += 1

        min_observed_distance = 0

        # In quick mode, stop early unless explicitly forced to evaluate all pairs.
        if not force_hamming_all:
            dists = [m["hamming_distance"] for m in matches]
            confs = [m["confidence"] for m in matches]
            return {
                "cell_barcode_1": cell1,
                "cell_barcode_2": cell2,
                "source_condition_1": source_1,
                "source_condition_2": source_2,
                "is_clone": len(matches) > 0,
                "num_matches": len(matches),
                "num_exact_matches": len(matches),
                "num_hamming_matches": 0,
                "num_barcode_comparisons": num_barcode_comparisons,
                "min_observed_distance": min_observed_distance,
                "best_hamming_distance": min(dists) if dists else np.nan,
                "mean_match_confidence": float(np.mean(confs)) if confs else np.nan,
                "matched_barcode_pairs": _serialize_matches(matches),
            }

    # Full pairwise barcode comparisons
    for e1 in entries1:
        for e2 in entries2:
            dist = calculate_hamming_distance(e1["barcode"], e2["barcode"])
            num_barcode_comparisons += 1
            if np.isnan(min_observed_distance) or dist < min_observed_distance:
                min_observed_distance = int(dist)

            if dist <= hamming_threshold:
                conf = calculate_match_confidence(
                    e1["umi_count"], e2["umi_count"], e1["rank"], e2["rank"]
                )
                matches.append(
                    {
                        "barcode_1": e1["barcode"],
                        "barcode_2": e2["barcode"],
                        "hamming_distance": int(dist),
                        "confidence": conf,
                        "umi_1": e1["umi_count"],
                        "umi_2": e2["umi_count"],
                        "rank_1": e1["rank"],
                        "rank_2": e2["rank"],
                    }
                )

    dists = [m["hamming_distance"] for m in matches]
    confs = [m["confidence"] for m in matches]
    num_exact = sum(1 for m in matches if m["hamming_distance"] == 0)
    num_hamming = sum(1 for m in matches if m["hamming_distance"] > 0)

    return {
        "cell_barcode_1": cell1,
        "cell_barcode_2": cell2,
        "source_condition_1": source_1,
        "source_condition_2": source_2,
        "is_clone": len(matches) > 0,
        "num_matches": len(matches),
        "num_exact_matches": num_exact,
        "num_hamming_matches": num_hamming,
        "num_barcode_comparisons": num_barcode_comparisons,
        "min_observed_distance": min_observed_distance,
        "best_hamming_distance": min(dists) if dists else np.nan,
        "mean_match_confidence": float(np.mean(confs)) if confs else np.nan,
        "matched_barcode_pairs": _serialize_matches(matches),
    }


def _init_within_worker(cell_data, cell_sets, params):
    """Initialize global state for within-condition workers."""
    global _WITHIN_CELL_DATA, _WITHIN_PARAMS
    _WITHIN_CELL_DATA = (cell_data, cell_sets)
    _WITHIN_PARAMS = params


def _within_worker(pair):
    """Worker for one within-condition cell pair."""
    cell1, cell2 = pair
    cell_data, cell_sets = _WITHIN_CELL_DATA
    return _compare_two_cells(
        cell1=cell1,
        cell2=cell2,
        entries1=cell_data[cell1],
        entries2=cell_data[cell2],
        set1=cell_sets[cell1],
        set2=cell_sets[cell2],
        hamming_threshold=_WITHIN_PARAMS["hamming_threshold"],
        quick_mode=_WITHIN_PARAMS["quick_mode"],
        force_hamming_all=_WITHIN_PARAMS["force_hamming_all"],
        source_1=_WITHIN_PARAMS["source_condition"],
        source_2=_WITHIN_PARAMS["source_condition"],
    )


def _init_cross_worker(cell_data_1, cell_sets_1, cell_data_2, cell_sets_2, params):
    """Initialize global state for cross-condition workers."""
    global _CROSS_CELL_DATA_1, _CROSS_CELL_DATA_2, _CROSS_PARAMS
    _CROSS_CELL_DATA_1 = (cell_data_1, cell_sets_1)
    _CROSS_CELL_DATA_2 = (cell_data_2, cell_sets_2)
    _CROSS_PARAMS = params


def _cross_worker(pair):
    """Worker for one cross-condition cell pair."""
    cell1, cell2 = pair
    cell_data_1, cell_sets_1 = _CROSS_CELL_DATA_1
    cell_data_2, cell_sets_2 = _CROSS_CELL_DATA_2
    return _compare_two_cells(
        cell1=cell1,
        cell2=cell2,
        entries1=cell_data_1[cell1],
        entries2=cell_data_2[cell2],
        set1=cell_sets_1[cell1],
        set2=cell_sets_2[cell2],
        hamming_threshold=_CROSS_PARAMS["hamming_threshold"],
        quick_mode=_CROSS_PARAMS["quick_mode"],
        force_hamming_all=_CROSS_PARAMS["force_hamming_all"],
        source_1=_CROSS_PARAMS["source_condition_1"],
        source_2=_CROSS_PARAMS["source_condition_2"],
    )


def _resolve_n_jobs(n_jobs):
    """Resolve worker count from user input."""
    if n_jobs == 0:
        return max(1, mp.cpu_count() - 1)
    return max(1, n_jobs)


def assign_clones_hamming(input_df, hamming_threshold=2, quick_mode=True, force_hamming_all=False, max_lineage_count=None, n_jobs=1, source_condition="input_1"):
    """
    Assign clones within one condition by pairwise cell comparisons.

    Parameters
    ----------
    input_df : pd.DataFrame
        Filtered lineage dataframe.
    hamming_threshold : int, optional
        Max Hamming distance for barcode matching.
    quick_mode : bool, optional
        Use exact-set intersection optimization.
    force_hamming_all : bool, optional
        If True, evaluate full barcode pairs even if exact overlap exists.
    max_lineage_count : int, optional
        Exclude cells with lineage_count >= this value.
    n_jobs : int, optional
        Number of worker processes (0 = auto).
    source_condition : str, optional
        Condition label used in detailed outputs.

    Returns
    -------
    pd.DataFrame
        Detailed pairwise comparison results with diagnostic columns including
        num_barcode_comparisons and min_observed_distance.
    """
    cell_data, cell_sets = _prepare_cell_barcode_data(input_df, max_lineage_count=max_lineage_count)
    cell_ids = sorted(cell_data.keys())

    if len(cell_ids) == 0:
        return pd.DataFrame(columns=DETAILED_PAIRWISE_COLUMNS)

    if len(cell_ids) == 1:
        cell = cell_ids[0]
        return pd.DataFrame(
            [
                {
                    "cell_barcode_1": cell,
                    "cell_barcode_2": cell,
                    "source_condition_1": source_condition,
                    "source_condition_2": source_condition,
                    "is_clone": False,
                    "num_matches": 0,
                    "num_exact_matches": 0,
                    "num_hamming_matches": 0,
                    "num_barcode_comparisons": 0,
                    "min_observed_distance": np.nan,
                    "best_hamming_distance": np.nan,
                    "mean_match_confidence": np.nan,
                    "matched_barcode_pairs": "",
                }
            ],
            columns=DETAILED_PAIRWISE_COLUMNS,
        )

    pair_iter = itertools.combinations(cell_ids, 2)
    total_pairs = len(cell_ids) * (len(cell_ids) - 1) // 2

    jobs = _resolve_n_jobs(n_jobs)
    params = {
        "hamming_threshold": hamming_threshold,
        "quick_mode": quick_mode,
        "force_hamming_all": force_hamming_all,
        "source_condition": source_condition,
    }

    print(f"Comparing {len(cell_ids)} cells ({total_pairs} pairs) with {jobs} worker(s) ...")

    if jobs > 1:
        with mp.Pool(
            jobs,
            initializer=_init_within_worker,
            initargs=(cell_data, cell_sets, params),
        ) as pool:
            rows = list(pool.imap(_within_worker, pair_iter, chunksize=200))
    else:
        _init_within_worker(cell_data, cell_sets, params)
        rows = [_within_worker(pair) for pair in pair_iter]

    return pd.DataFrame(rows, columns=DETAILED_PAIRWISE_COLUMNS)


def assign_clones_cross_condition(df1, df2, hamming_threshold=2, quick_mode=True, force_hamming_all=False, max_lineage_count=None, n_jobs=1, source_condition_1="input_1", source_condition_2="input_2"):
    """
    Assign clones across two conditions by all-against-all cell comparisons.

    Parameters
    ----------
    df1 : pd.DataFrame
        Filtered lineage dataframe for condition 1.
    df2 : pd.DataFrame
        Filtered lineage dataframe for condition 2.
    hamming_threshold : int, optional
        Max Hamming distance for barcode matching.
    quick_mode : bool, optional
        Use exact-set intersection optimization.
    force_hamming_all : bool, optional
        If True, evaluate full barcode pairs even if exact overlap exists.
    max_lineage_count : int, optional
        Exclude cells with lineage_count >= this value.
    n_jobs : int, optional
        Number of worker processes (0 = auto).
    source_condition_1 : str, optional
        Label for condition 1.
    source_condition_2 : str, optional
        Label for condition 2.

    Returns
    -------
    pd.DataFrame
        Detailed cross-condition pairwise comparison results with diagnostic
        columns including num_barcode_comparisons and min_observed_distance.
    """
    cell_data_1, cell_sets_1 = _prepare_cell_barcode_data(df1, max_lineage_count=max_lineage_count)
    cell_data_2, cell_sets_2 = _prepare_cell_barcode_data(df2, max_lineage_count=max_lineage_count)

    cells_1 = sorted(cell_data_1.keys())
    cells_2 = sorted(cell_data_2.keys())

    if len(cells_1) == 0 or len(cells_2) == 0:
        return pd.DataFrame(columns=DETAILED_PAIRWISE_COLUMNS)

    pair_iter = itertools.product(cells_1, cells_2)
    total_pairs = len(cells_1) * len(cells_2)

    jobs = _resolve_n_jobs(n_jobs)
    params = {
        "hamming_threshold": hamming_threshold,
        "quick_mode": quick_mode,
        "force_hamming_all": force_hamming_all,
        "source_condition_1": source_condition_1,
        "source_condition_2": source_condition_2,
    }

    print(
        f"Cross-condition comparison: {len(cells_1)} x {len(cells_2)} "
        f"cells ({total_pairs} pairs) with {jobs} worker(s) ..."
    )

    if jobs > 1:
        with mp.Pool(
            jobs,
            initializer=_init_cross_worker,
            initargs=(cell_data_1, cell_sets_1, cell_data_2, cell_sets_2, params),
        ) as pool:
            rows = list(pool.imap(_cross_worker, pair_iter, chunksize=200))
    else:
        _init_cross_worker(cell_data_1, cell_sets_1, cell_data_2, cell_sets_2, params)
        rows = [_cross_worker(pair) for pair in pair_iter]

    return pd.DataFrame(rows, columns=DETAILED_PAIRWISE_COLUMNS)


def generate_clone_groups(detailed_results):
    """
    Convert pairwise clone matches to clone groups via connected components.

    Parameters
    ----------
    detailed_results : pd.DataFrame
        Pairwise comparison dataframe generated by assign_clones_hamming/
        assign_clones_cross_condition.

    Returns
    -------
    pd.DataFrame
        Dataframe with cloneID assignments per cell.
    """
    if detailed_results.empty:
        return pd.DataFrame(columns=["cloneID", "source_condition", "cell_barcode", "clone_size"])

    # Build graph over (source_condition, cell_barcode) nodes.
    all_nodes = set()
    adjacency = defaultdict(set)

    for row in detailed_results.itertuples(index=False):
        node1 = (row.source_condition_1, row.cell_barcode_1)
        node2 = (row.source_condition_2, row.cell_barcode_2)
        all_nodes.add(node1)
        all_nodes.add(node2)

        if bool(row.is_clone):
            adjacency[node1].add(node2)
            adjacency[node2].add(node1)

    visited = set()
    components = []

    for node in sorted(all_nodes):
        if node in visited:
            continue
        stack = [node]
        component = []
        visited.add(node)

        while stack:
            cur = stack.pop()
            component.append(cur)
            for nei in adjacency[cur]:
                if nei not in visited:
                    visited.add(nei)
                    stack.append(nei)

        components.append(sorted(component))

    components = sorted(components, key=lambda x: len(x), reverse=True)

    rows = []
    for idx, comp in enumerate(components, start=1):
        for source_condition, cell_barcode in comp:
            rows.append(
                {
                    "cloneID": idx,
                    "source_condition": source_condition,
                    "cell_barcode": cell_barcode,
                    "clone_size": len(comp),
                }
            )

    return pd.DataFrame(rows)


def extract_unique_matched_cells(detailed_results, source_condition):
    """
    Extract unique matched cell barcodes for a given source condition.

    Parameters
    ----------
    detailed_results : pd.DataFrame
        Pairwise comparison dataframe.
    source_condition : str
        Source condition label to extract matched cells for.

    Returns
    -------
    pd.DataFrame
        Single-column dataframe of unique matched cell barcodes.
    """
    if detailed_results.empty:
        return pd.DataFrame(columns=["cell_barcode"])

    clone_rows = detailed_results[detailed_results["is_clone"]].copy()
    if clone_rows.empty:
        return pd.DataFrame(columns=["cell_barcode"])

    cells_side_1 = clone_rows.loc[
        clone_rows["source_condition_1"] == source_condition, "cell_barcode_1"
    ]
    cells_side_2 = clone_rows.loc[
        clone_rows["source_condition_2"] == source_condition, "cell_barcode_2"
    ]

    unique_cells = sorted(set(cells_side_1).union(set(cells_side_2)))
    return pd.DataFrame({"cell_barcode": unique_cells})


def _build_clone_barcode_outputs(filtered_df, clone_groups_df):
    """Build legacy cloneID output tables from clone groups."""
    if clone_groups_df.empty:
        empty_clonebarcode = pd.DataFrame(columns=["cloneID", "clone_barcode"])
        empty_cellbarcode = pd.DataFrame(columns=["cloneID", "cell_barcode"])
        empty_counts = pd.DataFrame(columns=["clone_barcode", "count", "cloneID"])
        return empty_clonebarcode, empty_cellbarcode, empty_counts

    # Map each cell to unique lineage barcodes from filtered dataframe.
    cell_to_barcodes = (
        filtered_df.groupby("cell_barcode")["corrected_lineage_barcode"]
        .apply(lambda x: sorted(set(x.astype(str))))
        .to_dict()
    )

    clone_barcodes = []
    cloneid_cellbarcode_rows = []

    for clone_id, grp in clone_groups_df.groupby("cloneID"):
        cells = sorted(grp["cell_barcode"].tolist())
        union_barcodes = sorted(
            {
                bc
                for cell in cells
                for bc in cell_to_barcodes.get(cell, [])
            }
        )

        clone_barcode = "_".join(union_barcodes) if union_barcodes else ""
        clone_barcodes.append({"cloneID": int(clone_id), "clone_barcode": clone_barcode})

        for cell in cells:
            cloneid_cellbarcode_rows.append({"cloneID": int(clone_id), "cell_barcode": cell})

    cloneid_clonebarcode_df = pd.DataFrame(clone_barcodes).sort_values("cloneID")
    cloneid_cellbarcode_df = pd.DataFrame(cloneid_cellbarcode_rows).sort_values(
        ["cloneID", "cell_barcode"]
    )

    clone_counts = (
        cloneid_cellbarcode_df.groupby("cloneID")
        .size()
        .reset_index(name="count")
        .merge(cloneid_clonebarcode_df, on="cloneID", how="left")
        [["clone_barcode", "count", "cloneID"]]
        .sort_values(["count", "cloneID"], ascending=[False, True])
        .reset_index(drop=True)
    )

    return cloneid_clonebarcode_df, cloneid_cellbarcode_df, clone_counts


def save_results(
    filebase,
    z_threshold,
    max_lineage_count,
    df,
    filtered_df,
    cloneid_clonebarcode_df,
    cloneid_cellbarcode_df,
    detailed_results_df,
    clone_groups_df,
    matched_cells_df,
    qc_results,
):
    """Save analysis outputs, including new clone assignment artifacts."""
    if max_lineage_count:
        zcfilebase = f"{filebase}.z{z_threshold}_mlc{max_lineage_count}"
        zfilebase = f"{filebase}.z{z_threshold}"
    else:
        zfilebase = f"{filebase}.z{z_threshold}"
        zcfilebase = f"{filebase}.z{z_threshold}"

    df.to_csv(f"{filebase}.rewind_statistics.tsv", sep="\t", index=False)
    filtered_df.to_csv(f"{zfilebase}.rewind_filtered.tsv", sep="\t", index=False)
    cloneid_clonebarcode_df.to_csv(
        f"{zcfilebase}.cloneID_cloneBarcode.tsv", sep="\t", index=False
    )
    cloneid_cellbarcode_df.to_csv(
        f"{zcfilebase}.cloneID_cellBarcode.tsv", sep="\t", index=False
    )

    # New outputs
    detailed_results_for_write = detailed_results_df.copy()
    for col in DETAILED_PAIRWISE_COLUMNS:
        if col not in detailed_results_for_write.columns:
            detailed_results_for_write[col] = np.nan
    detailed_results_for_write = detailed_results_for_write[DETAILED_PAIRWISE_COLUMNS]

    detailed_results_for_write.to_csv(
        f"{zcfilebase}.clone_pairwise_detailed.tsv", sep="\t", index=False
    )
    clone_groups_df.to_csv(f"{zcfilebase}.clone_groups.tsv", sep="\t", index=False)
    matched_cells_df.to_csv(f"{zcfilebase}.matched_cells.tsv", sep="\t", index=False)

    with open(f"{zfilebase}.truncated_Poisson_fit_results.txt", "w") as f:
        for key, value in qc_results.items():
            f.write(f"{key}: {value}\n")


def plot_lineage_barcode_qc(filebase, empirical_counts, expected_counts, expected_cdf, empirical_cdf):
    """
    Generate QC plots for lineage barcode count distribution fitting.

    Parameters
    ----------
    filebase : str
        File prefix for output plot.
    empirical_counts : np.ndarray
        Empirical counts per lineage-count bucket.
    expected_counts : np.ndarray
        Expected counts from zero-truncated Poisson fit.
    expected_cdf : np.ndarray
        Expected CDF values.
    empirical_cdf : np.ndarray
        Empirical CDF values.
    """
    x_values = np.arange(1, len(expected_counts) + 1)

    plot_data = {
        1: {
            "y_expected": expected_counts,
            "y_empirical": empirical_counts,
            "ylabel": "Counts",
            "title": "Expected vs Empirical Counts",
        },
        2: {
            "y_expected": expected_counts / expected_counts.sum(),
            "y_empirical": empirical_counts / empirical_counts.sum(),
            "ylabel": "Fraction",
            "title": "Expected vs Empirical Fractions",
        },
        3: {
            "y_expected": expected_cdf,
            "y_empirical": empirical_cdf,
            "ylabel": "CDF",
            "title": "Expected vs Empirical CDFs",
        },
    }

    plt.figure(figsize=(18, 6))

    for i in range(1, 4):
        plt.subplot(1, 3, i)
        plt.plot(x_values, plot_data[i]["y_expected"], label="Expected")
        plt.plot(x_values, plot_data[i]["y_empirical"], label="Empirical")
        plt.xlabel("Number of Lineages per Cell")
        plt.ylabel(plot_data[i]["ylabel"])
        plt.legend()
        plt.title(plot_data[i]["title"])

    plt.tight_layout()
    plt.savefig(f"{filebase}.truncated_Poisson_fit_plots.png")
    plt.close()


def plot_clone_barcode_qc(filebase, clone_counts):
    """
    Plot distribution of number of cells per clone.

    Parameters
    ----------
    filebase : str
        File prefix for output plot.
    clone_counts : pd.DataFrame
        DataFrame containing at least the 'count' column.
    """
    if clone_counts.empty:
        print("No clone counts available for clone barcode QC plot; skipping.")
        return

    data = clone_counts["count"]
    plt.figure(figsize=(6, 6))
    sns.histplot(data, bins=range(1, int(max(data)) + 2), discrete=True, stat="density")

    ax = plt.gca()
    ax.set_xticks(range(1, int(max(data)) + 1))
    ax.set_xticklabels(range(1, int(max(data)) + 1))
    ax.set_xlabel("Number of Cells per Clone")
    ax.set_ylabel("Fraction of Clones")

    plt.tight_layout()
    plt.savefig(f"{filebase}.clone_barcode_distribution_qc_plot.png")
    plt.close()


def main():
    """Entry point: filter noise, assign clones, and write outputs."""
    args = parse_args()

    if args.cross_condition and not args.input_file_2:
        raise ValueError("--cross-condition requires --input-file-2")

    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        prefix = os.path.basename(args.input_file).rsplit(".", 1)[0]
        filebase = os.path.join(args.output_dir, prefix)
    else:
        filebase = args.input_file.rsplit(".", 1)[0]

    if args.max_lineage_count:
        zcfilebase = f"{filebase}.z{args.z_threshold}_mlc{args.max_lineage_count}"
        zfilebase = f"{filebase}.z{args.z_threshold}"
    else:
        zfilebase = f"{filebase}.z{args.z_threshold}"
        zcfilebase = f"{filebase}.z{args.z_threshold}"

    print(filebase)
    print(zfilebase)
    print(zcfilebase)

    # Load primary input and run noise filtering/QC
    df = pd.read_csv(args.input_file, sep="\t")
    df, filtered_df = filter_noise(df, args.top_n, args.z_threshold, args.min_group_size)

    lineage_count_df = filtered_df[["cell_barcode", "lineage_count"]].drop_duplicates().reset_index(drop=True)
    empirical_counts, expected_counts, empirical_cdf, expected_cdf, results = test_truncated_poisson_distribution(
        lineage_count_df, "lineage_count"
    )
    plot_lineage_barcode_qc(zfilebase, empirical_counts, expected_counts, expected_cdf, empirical_cdf)

    # Clone assignment mode
    if args.cross_condition:
        df2 = pd.read_csv(args.input_file_2, sep="\t")
        _, filtered_df_2 = filter_noise(df2, args.top_n, args.z_threshold, args.min_group_size)

        detailed_results_df = assign_clones_cross_condition(
            filtered_df,
            filtered_df_2,
            hamming_threshold=args.hamming_threshold,
            quick_mode=args.quick_mode,
            force_hamming_all=args.force_hamming_all,
            max_lineage_count=args.max_lineage_count,
            n_jobs=args.n_jobs,
            source_condition_1=os.path.basename(args.input_file),
            source_condition_2=os.path.basename(args.input_file_2),
        )

        clone_groups_df = generate_clone_groups(detailed_results_df)

        matched_1 = extract_unique_matched_cells(detailed_results_df, os.path.basename(args.input_file))
        matched_1["source_condition"] = os.path.basename(args.input_file)

        matched_2 = extract_unique_matched_cells(detailed_results_df, os.path.basename(args.input_file_2))
        matched_2["source_condition"] = os.path.basename(args.input_file_2)

        matched_cells_df = pd.concat([matched_1, matched_2], ignore_index=True)
        matched_cells_df = matched_cells_df[["source_condition", "cell_barcode"]]

        # Keep legacy outputs tied to input-file cells for compatibility.
        cloneid_clonebarcode_df, cloneid_cellbarcode_df, clone_counts = _build_clone_barcode_outputs(
            filtered_df, clone_groups_df[clone_groups_df["source_condition"] == os.path.basename(args.input_file)]
        )
    else:
        detailed_results_df = assign_clones_hamming(
            filtered_df,
            hamming_threshold=args.hamming_threshold,
            quick_mode=args.quick_mode,
            force_hamming_all=args.force_hamming_all,
            max_lineage_count=args.max_lineage_count,
            n_jobs=args.n_jobs,
            source_condition=os.path.basename(args.input_file),
        )

        clone_groups_df = generate_clone_groups(detailed_results_df)
        matched_cells_df = extract_unique_matched_cells(
            detailed_results_df, os.path.basename(args.input_file)
        )

        cloneid_clonebarcode_df, cloneid_cellbarcode_df, clone_counts = _build_clone_barcode_outputs(
            filtered_df, clone_groups_df
        )

    plot_clone_barcode_qc(zcfilebase, clone_counts)

    save_results(
        filebase=filebase,
        z_threshold=args.z_threshold,
        max_lineage_count=args.max_lineage_count,
        df=df,
        filtered_df=filtered_df,
        cloneid_clonebarcode_df=cloneid_clonebarcode_df,
        cloneid_cellbarcode_df=cloneid_cellbarcode_df,
        detailed_results_df=detailed_results_df,
        clone_groups_df=clone_groups_df,
        matched_cells_df=matched_cells_df,
        qc_results=results,
    )


if __name__ == "__main__":
    main()
