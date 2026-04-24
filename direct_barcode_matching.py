#!/usr/bin/env python3
"""
Direct Cell-to-Cell Barcode Matching for RewindPipeline.

This script matches T1 cells (Ctrl and/or Treatment) to T0 cells directly from
Step 9 ``*_corrected_UMI_counts.tsv`` outputs, bypassing clone assignment.

Input file format (TSV, header required):
    - cell_barcode
    - corrected_lineage_barcode
    - umi_count
    - total_molecules_per_cell
    - total_lineages_per_cell

Output:
    One CSV per T1 condition in the output directory.

Example:
    python direct_barcode_matching.py \
      --t0 T0_corrected_UMI_counts.tsv \
      --t1-ctrl T1_Ctrl_corrected_UMI_counts.tsv \
      --t1-treatment T1_Treatment_corrected_UMI_counts.tsv \
      --output-dir results/ \
      --hamming-threshold 2 \
      --jaccard-threshold 0.6 \
      --min-matched-barcodes 2 \
      --use-umi-weighting \
      --n-jobs 4
"""

from __future__ import annotations

import argparse
import csv
import logging
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

REQUIRED_COLUMNS = {
    "cell_barcode",
    "corrected_lineage_barcode",
    "umi_count",
    "total_molecules_per_cell",
    "total_lineages_per_cell",
}

NO_MATCH_PLACEHOLDER = "NO_MATCH"


@dataclass
class CellBarcodeData:
    """Container for per-cell barcode and UMI information."""

    cell_to_barcodes: Dict[str, Set[str]]
    cell_to_barcode_umi: Dict[str, Dict[str, int]]
    cell_total_umis: Dict[str, int]
    cell_total_lineages: Dict[str, int]


@dataclass
class MatchResult:
    """Stores one matched candidate between a T1 cell and a T0 cell."""

    t0_cell_barcode: str
    t0_barcodes: str
    t0_barcode_count: int
    t1_cell_barcode: str
    t1_barcodes: str
    t1_barcode_count: int
    matched_barcodes: int
    jaccard_similarity: float
    umi_weighted_similarity: Optional[float]
    match_rank: int


# Globals used by multiprocessing workers.
WORKER_REF = {}


def setup_logging(level: str) -> None:
    """Configure process-wide logging."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def positive_int(value: str) -> int:
    """Argparse validator for strictly positive integers."""
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"Expected a positive integer, got {value!r}")
    return parsed


def non_negative_int(value: str) -> int:
    """Argparse validator for non-negative integers."""
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"Expected a non-negative integer, got {value!r}")
    return parsed


def bounded_float_0_1(value: str) -> float:
    """Argparse validator for floats in [0, 1]."""
    parsed = float(value)
    if parsed < 0.0 or parsed > 1.0:
        raise argparse.ArgumentTypeError(f"Expected a value in [0, 1], got {value!r}")
    return parsed


def calculate_hamming_distance(barcode_a: str, barcode_b: str) -> int:
    """Return Hamming distance between two equal-length barcodes.

    Raises
    ------
    ValueError
        If barcode lengths are not equal.
    """
    if len(barcode_a) != len(barcode_b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(ch1 != ch2 for ch1, ch2 in zip(barcode_a, barcode_b))


def parse_step9_tsv(path: Path) -> CellBarcodeData:
    """Parse Step 9 TSV into cell-centric dictionaries with validation."""
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if not path.is_file():
        raise ValueError(f"Input path is not a file: {path}")

    cell_to_barcodes: Dict[str, Set[str]] = defaultdict(set)
    cell_to_barcode_umi: Dict[str, Dict[str, int]] = defaultdict(dict)
    cell_total_umis: Dict[str, int] = {}
    cell_total_lineages: Dict[str, int] = {}

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header in TSV file: {path}")

        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"TSV missing required columns in {path}: {sorted(missing)}"
            )

        row_count = 0
        for row_idx, row in enumerate(reader, start=2):
            row_count += 1
            cell = (row.get("cell_barcode") or "").strip()
            barcode = (row.get("corrected_lineage_barcode") or "").strip()
            umi_raw = (row.get("umi_count") or "").strip()
            total_umi_raw = (row.get("total_molecules_per_cell") or "").strip()
            total_lin_raw = (row.get("total_lineages_per_cell") or "").strip()

            if not cell:
                raise ValueError(f"Empty cell_barcode at {path}:{row_idx}")
            if not barcode:
                raise ValueError(f"Empty corrected_lineage_barcode at {path}:{row_idx}")

            try:
                umi_count = int(umi_raw)
                total_umi = int(total_umi_raw)
                total_lineages = int(total_lin_raw)
            except ValueError as exc:
                raise ValueError(
                    f"Non-integer count at {path}:{row_idx} | row={row}"
                ) from exc

            if umi_count < 0 or total_umi < 0 or total_lineages < 0:
                raise ValueError(f"Negative counts are not allowed at {path}:{row_idx}")

            if barcode in cell_to_barcode_umi[cell]:
                cell_to_barcode_umi[cell][barcode] += umi_count
            else:
                cell_to_barcode_umi[cell][barcode] = umi_count
            cell_to_barcodes[cell].add(barcode)

            if cell in cell_total_umis and cell_total_umis[cell] != total_umi:
                logging.warning(
                    "Inconsistent total_molecules_per_cell for cell %s in %s; "
                    "using sum from barcode rows.",
                    cell,
                    path,
                )
            if cell in cell_total_lineages and cell_total_lineages[cell] != total_lineages:
                logging.warning(
                    "Inconsistent total_lineages_per_cell for cell %s in %s; "
                    "using observed unique barcodes.",
                    cell,
                    path,
                )

            cell_total_umis[cell] = total_umi
            cell_total_lineages[cell] = total_lineages

    # Normalize totals from observed barcode UMIs to guarantee consistency.
    for cell, barcode_to_umi in cell_to_barcode_umi.items():
        observed_total = sum(barcode_to_umi.values())
        observed_lineages = len(cell_to_barcodes[cell])
        cell_total_umis[cell] = observed_total
        cell_total_lineages[cell] = observed_lineages

    if row_count == 0:
        logging.warning("Input file is empty (header only): %s", path)

    return CellBarcodeData(
        cell_to_barcodes=dict(cell_to_barcodes),
        cell_to_barcode_umi={k: dict(v) for k, v in cell_to_barcode_umi.items()},
        cell_total_umis=cell_total_umis,
        cell_total_lineages=cell_total_lineages,
    )


def build_reference_index(
    t0_data: CellBarcodeData,
) -> Tuple[Dict[int, List[str]], Dict[str, Set[str]], Dict[str, int]]:
    """Build index structures for efficient barcode-to-T0 candidate lookup.

    Returns
    -------
    barcodes_by_length
        Mapping of barcode length to all unique T0 barcodes with that length.
    barcode_to_t0_cells
        Mapping barcode -> set of T0 cell barcodes containing it.
    t0_barcode_umi_total
        Barcode-level total UMI across all T0 cells (optional diagnostics/weighting).
    """
    barcodes_by_length: Dict[int, List[str]] = defaultdict(list)
    barcode_to_t0_cells: Dict[str, Set[str]] = defaultdict(set)
    t0_barcode_umi_total: Dict[str, int] = defaultdict(int)

    seen: Set[str] = set()
    for t0_cell, barcodes in t0_data.cell_to_barcodes.items():
        barcode_umi_map = t0_data.cell_to_barcode_umi.get(t0_cell, {})
        for barcode in barcodes:
            barcode_to_t0_cells[barcode].add(t0_cell)
            t0_barcode_umi_total[barcode] += barcode_umi_map.get(barcode, 0)
            if barcode not in seen:
                barcodes_by_length[len(barcode)].append(barcode)
                seen.add(barcode)

    return dict(barcodes_by_length), dict(barcode_to_t0_cells), dict(t0_barcode_umi_total)


def perform_fuzzy_barcode_matching(
    query_barcodes: Iterable[str],
    ref_barcodes_by_length: Dict[int, List[str]],
    hamming_threshold: int,
) -> Dict[str, List[Tuple[str, int]]]:
    """Match each query barcode to reference barcodes within Hamming threshold.

    Only barcodes with equal length are compared (Hamming-compatible).
    Returns mapping:
        query_barcode -> [(ref_barcode, hamming_distance), ...]
    """
    results: Dict[str, List[Tuple[str, int]]] = {}

    for query in query_barcodes:
        refs = ref_barcodes_by_length.get(len(query), [])
        matches: List[Tuple[str, int]] = []

        for ref in refs:
            # Inline distance with early stop for speed on large datasets.
            mismatches = 0
            for ch1, ch2 in zip(query, ref):
                if ch1 != ch2:
                    mismatches += 1
                    if mismatches > hamming_threshold:
                        break
            if mismatches <= hamming_threshold:
                matches.append((ref, mismatches))

        results[query] = matches

    return results


def select_one_to_one_barcode_matches(
    barcode_matches: Dict[str, List[Tuple[str, int]]],
    allowed_t0_barcodes: Set[str],
) -> Tuple[List[Tuple[str, str]], List[int]]:
    """Select one-to-one T1<->T0 barcode matches using greedy bipartite matching.

    We first materialize all candidate edges (t1_barcode, t0_barcode, distance), then
    sort by ascending Hamming distance and deterministic lexical tiebreakers.
    The greedy walk keeps an edge only when both barcodes are still unmatched.
    This prevents many-to-one inflation where one barcode could otherwise match
    multiple partners and incorrectly increase intersection/Jaccard values.
    """
    candidate_edges: List[Tuple[int, str, str]] = []
    for t1_bc, matched_refs in barcode_matches.items():
        for t0_bc, distance in matched_refs:
            if t0_bc in allowed_t0_barcodes:
                candidate_edges.append((distance, t1_bc, t0_bc))

    candidate_edges.sort(key=lambda item: (item[0], item[1], item[2]))

    used_t1: Set[str] = set()
    used_t0: Set[str] = set()
    selected_pairs: List[Tuple[str, str]] = []
    selected_distances: List[int] = []

    for distance, t1_bc, t0_bc in candidate_edges:
        if t1_bc in used_t1 or t0_bc in used_t0:
            continue
        used_t1.add(t1_bc)
        used_t0.add(t0_bc)
        selected_pairs.append((t1_bc, t0_bc))
        selected_distances.append(distance)

    return selected_pairs, selected_distances


def calculate_jaccard_similarity(
    t0_barcodes: Set[str],
    t1_barcodes: Set[str],
    matched_pairs: Sequence[Tuple[str, str]],
    t0_umi_by_barcode: Optional[Dict[str, int]] = None,
    t1_umi_by_barcode: Optional[Dict[str, int]] = None,
) -> Tuple[float, Optional[float], int]:
    """Compute standard and optional UMI-weighted fuzzy Jaccard similarity.

    Parameters
    ----------
    matched_pairs
        Sequence of selected (t1_barcode, t0_barcode) matches.
    """
    matched_count = len(matched_pairs)
    union_size = len(t0_barcodes) + len(t1_barcodes) - matched_count

    if union_size <= 0:
        base_jaccard = 0.0
    else:
        base_jaccard = matched_count / union_size

    if t0_umi_by_barcode is None or t1_umi_by_barcode is None:
        return base_jaccard, None, union_size

    # Weighted intersection from matched barcode pairs.
    matched_weight = 0
    for t1_bc, t0_bc in matched_pairs:
        matched_weight += min(
            t1_umi_by_barcode.get(t1_bc, 0),
            t0_umi_by_barcode.get(t0_bc, 0),
        )

    t0_total = sum(t0_umi_by_barcode.values())
    t1_total = sum(t1_umi_by_barcode.values())
    weighted_union = t0_total + t1_total - matched_weight

    if weighted_union <= 0:
        weighted_jaccard = 0.0
    else:
        weighted_jaccard = matched_weight / weighted_union

    return base_jaccard, weighted_jaccard, union_size


def init_worker(worker_ref: dict) -> None:
    """Initialize multiprocessing worker with read-only matching references."""
    global WORKER_REF
    WORKER_REF = worker_ref


def match_single_t1_cell(
    t1_cell: str,
    t1_barcodes: Set[str],
    t1_umi_by_barcode: Dict[str, int],
    min_matched_barcodes: int,
    jaccard_threshold: float,
    use_umi_weighting: bool,
) -> List[MatchResult]:
    """Match one T1 cell against all candidate T0 cells and return ranked matches."""
    ref_barcodes_by_length = WORKER_REF["ref_barcodes_by_length"]
    barcode_to_t0_cells = WORKER_REF["barcode_to_t0_cells"]
    t0_cell_to_barcodes = WORKER_REF["t0_cell_to_barcodes"]
    t0_cell_to_umi = WORKER_REF["t0_cell_to_umi"]
    hamming_threshold = WORKER_REF["hamming_threshold"]

    barcode_matches = perform_fuzzy_barcode_matching(
        query_barcodes=t1_barcodes,
        ref_barcodes_by_length=ref_barcodes_by_length,
        hamming_threshold=hamming_threshold,
    )

    candidate_t0_cells: Set[str] = set()
    for matched_refs in barcode_matches.values():
        for t0_bc, _ in matched_refs:
            candidate_t0_cells.update(barcode_to_t0_cells.get(t0_bc, set()))

    if not candidate_t0_cells:
        return []

    scored_candidates: List[Tuple[float, int, float, MatchResult]] = []

    for t0_cell in candidate_t0_cells:
        t0_barcodes = t0_cell_to_barcodes[t0_cell]
        t0_umi_by_barcode = t0_cell_to_umi.get(t0_cell, {})

        selected_pairs, selected_distances = select_one_to_one_barcode_matches(
            barcode_matches=barcode_matches,
            allowed_t0_barcodes=t0_barcodes,
        )

        matched_count = len(selected_pairs)
        if matched_count < min_matched_barcodes:
            continue

        jaccard, umi_weighted, _ = calculate_jaccard_similarity(
            t0_barcodes=t0_barcodes,
            t1_barcodes=t1_barcodes,
            matched_pairs=selected_pairs,
            t0_umi_by_barcode=t0_umi_by_barcode if use_umi_weighting else None,
            t1_umi_by_barcode=t1_umi_by_barcode if use_umi_weighting else None,
        )

        score_for_threshold = umi_weighted if (use_umi_weighting and umi_weighted is not None) else jaccard
        if score_for_threshold < jaccard_threshold:
            continue

        mean_distance = (sum(selected_distances) / len(selected_distances)) if selected_distances else float("inf")
        result = MatchResult(
            t0_cell_barcode=t0_cell,
            t0_barcodes=";".join(sorted(t0_barcodes)),
            t0_barcode_count=len(t0_barcodes),
            t1_cell_barcode=t1_cell,
            t1_barcodes=";".join(sorted(t1_barcodes)),
            t1_barcode_count=len(t1_barcodes),
            matched_barcodes=matched_count,
            jaccard_similarity=round(jaccard, 6),
            umi_weighted_similarity=round(umi_weighted, 6) if umi_weighted is not None else None,
            match_rank=0,
        )

        scored_candidates.append((score_for_threshold, matched_count, mean_distance, result))

    # Ranking: score desc, matched count desc, mean distance asc, lexical T0 cell
    scored_candidates.sort(
        key=lambda item: (-item[0], -item[1], item[2], item[3].t0_cell_barcode)
    )

    ranked_results: List[MatchResult] = []
    for rank, (_, _, _, result) in enumerate(scored_candidates, start=1):
        result.match_rank = rank
        ranked_results.append(result)

    return ranked_results


def write_match_csv(
    output_path: Path,
    rows: Sequence[MatchResult],
    include_umi_weighted: bool,
) -> None:
    """Write match results to CSV with deterministic column order."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "T0_cell_barcode",
        "T0_barcodes",
        "T0_barcode_count",
        "T1_cell_barcode",
        "T1_barcodes",
        "T1_barcode_count",
        "matched_barcodes",
        "jaccard_similarity",
    ]
    if include_umi_weighted:
        fieldnames.append("umi_weighted_similarity")
    fieldnames.append("match_rank")

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            payload = {
                "T0_cell_barcode": row.t0_cell_barcode,
                "T0_barcodes": row.t0_barcodes,
                "T0_barcode_count": row.t0_barcode_count,
                "T1_cell_barcode": row.t1_cell_barcode,
                "T1_barcodes": row.t1_barcodes,
                "T1_barcode_count": row.t1_barcode_count,
                "matched_barcodes": row.matched_barcodes,
                "jaccard_similarity": f"{row.jaccard_similarity:.6f}",
                "match_rank": row.match_rank,
            }
            if include_umi_weighted:
                payload["umi_weighted_similarity"] = (
                    f"{row.umi_weighted_similarity:.6f}"
                    if row.umi_weighted_similarity is not None
                    else ""
                )
            writer.writerow(payload)


def log_summary_statistics(
    label: str,
    all_rows: Sequence[MatchResult],
    total_t1_cells: int,
) -> None:
    """Log summary statistics for one comparison run."""
    matched_t1_cells = len({r.t1_cell_barcode for r in all_rows if r.matched_barcodes > 0})
    unmatched_t1_cells = max(total_t1_cells - matched_t1_cells, 0)
    jaccard_scores = [
        r.jaccard_similarity
        for r in all_rows
        if r.match_rank == 1 and r.matched_barcodes > 0
    ]
    weighted_scores = [
        r.umi_weighted_similarity
        for r in all_rows
        if r.match_rank == 1 and r.matched_barcodes > 0 and r.umi_weighted_similarity is not None
    ]

    logging.info("--- Summary: %s ---", label)
    logging.info("Total T1 cells: %d", total_t1_cells)
    logging.info("Matched T1 cells: %d", matched_t1_cells)
    logging.info("Unmatched T1 cells: %d", unmatched_t1_cells)
    logging.info("Total emitted match rows: %d", len(all_rows))

    if jaccard_scores:
        logging.info(
            "Best-match Jaccard | mean=%.4f median=%.4f min=%.4f max=%.4f",
            mean(jaccard_scores),
            median(jaccard_scores),
            min(jaccard_scores),
            max(jaccard_scores),
        )
    else:
        logging.info("No passing matches, so no Jaccard summary available.")

    if weighted_scores:
        logging.info(
            "Best-match UMI-weighted | mean=%.4f median=%.4f min=%.4f max=%.4f",
            mean(weighted_scores),
            median(weighted_scores),
            min(weighted_scores),
            max(weighted_scores),
        )


def _run_condition_matching(
    condition_label: str,
    t0_data: CellBarcodeData,
    t1_data: CellBarcodeData,
    output_path: Path,
    hamming_threshold: int,
    jaccard_threshold: float,
    min_matched_barcodes: int,
    use_umi_weighting: bool,
    n_jobs: int,
    top_n: int,
) -> None:
    """Run matching for one T1 condition and write output CSV."""
    logging.info("\n=== Matching condition: %s ===", condition_label)
    logging.info("T0 cells: %d | T1 cells: %d", len(t0_data.cell_to_barcodes), len(t1_data.cell_to_barcodes))

    ref_barcodes_by_length, barcode_to_t0_cells, _ = build_reference_index(t0_data)
    total_ref_barcodes = sum(len(v) for v in ref_barcodes_by_length.values())
    logging.info("T0 unique barcodes indexed: %d", total_ref_barcodes)

    worker_ref = {
        "ref_barcodes_by_length": ref_barcodes_by_length,
        "barcode_to_t0_cells": barcode_to_t0_cells,
        "t0_cell_to_barcodes": t0_data.cell_to_barcodes,
        "t0_cell_to_umi": t0_data.cell_to_barcode_umi,
        "hamming_threshold": hamming_threshold,
    }

    t1_items = list(t1_data.cell_to_barcodes.items())
    total_t1 = len(t1_items)

    rows: List[MatchResult] = []
    started = time.time()

    if total_t1 == 0:
        logging.warning("No T1 cells found for condition %s; writing empty output CSV.", condition_label)
        write_match_csv(output_path, rows, include_umi_weighted=use_umi_weighting)
        return

    if n_jobs > 1:
        logging.info("Running matching with multiprocessing (%d workers)", n_jobs)
        with mp.Pool(processes=n_jobs, initializer=init_worker, initargs=(worker_ref,)) as pool:
            args_iter = (
                (
                    t1_cell,
                    t1_barcodes,
                    t1_data.cell_to_barcode_umi.get(t1_cell, {}),
                    min_matched_barcodes,
                    jaccard_threshold,
                    use_umi_weighting,
                )
                for t1_cell, t1_barcodes in t1_items
            )

            for idx, ranked_matches in enumerate(pool.starmap(match_single_t1_cell, args_iter, chunksize=50), start=1):
                t1_cell, t1_barcodes = t1_items[idx - 1]
                if ranked_matches:
                    if top_n > 0:
                        rows.extend(ranked_matches[:top_n])
                    else:
                        rows.extend(ranked_matches)
                else:
                    rows.append(
                        MatchResult(
                            t0_cell_barcode=NO_MATCH_PLACEHOLDER,
                            t0_barcodes=NO_MATCH_PLACEHOLDER,
                            t0_barcode_count=0,
                            t1_cell_barcode=t1_cell,
                            t1_barcodes=";".join(sorted(t1_barcodes)),
                            t1_barcode_count=len(t1_barcodes),
                            matched_barcodes=0,
                            jaccard_similarity=0.0,
                            umi_weighted_similarity=0.0 if use_umi_weighting else None,
                            match_rank=1,
                        )
                    )
                if idx % 200 == 0 or idx == total_t1:
                    elapsed = time.time() - started
                    logging.info(
                        "Progress: %d/%d T1 cells (%.1f%%) | elapsed %.1fs",
                        idx,
                        total_t1,
                        100.0 * idx / total_t1,
                        elapsed,
                    )
    else:
        logging.info("Running matching in single-process mode")
        init_worker(worker_ref)
        for idx, (t1_cell, t1_barcodes) in enumerate(t1_items, start=1):
            ranked_matches = match_single_t1_cell(
                t1_cell=t1_cell,
                t1_barcodes=t1_barcodes,
                t1_umi_by_barcode=t1_data.cell_to_barcode_umi.get(t1_cell, {}),
                min_matched_barcodes=min_matched_barcodes,
                jaccard_threshold=jaccard_threshold,
                use_umi_weighting=use_umi_weighting,
            )
            if ranked_matches:
                if top_n > 0:
                    rows.extend(ranked_matches[:top_n])
                else:
                    rows.extend(ranked_matches)
            else:
                rows.append(
                    MatchResult(
                        t0_cell_barcode=NO_MATCH_PLACEHOLDER,
                        t0_barcodes=NO_MATCH_PLACEHOLDER,
                        t0_barcode_count=0,
                        t1_cell_barcode=t1_cell,
                        t1_barcodes=";".join(sorted(t1_barcodes)),
                        t1_barcode_count=len(t1_barcodes),
                        matched_barcodes=0,
                        jaccard_similarity=0.0,
                        umi_weighted_similarity=0.0 if use_umi_weighting else None,
                        match_rank=1,
                    )
                )
            if idx % 200 == 0 or idx == total_t1:
                elapsed = time.time() - started
                logging.info(
                    "Progress: %d/%d T1 cells (%.1f%%) | elapsed %.1fs",
                    idx,
                    total_t1,
                    100.0 * idx / total_t1,
                    elapsed,
                )

    write_match_csv(output_path, rows, include_umi_weighted=use_umi_weighting)
    log_summary_statistics(condition_label, rows, total_t1_cells=total_t1)
    logging.info("Output written: %s", output_path)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Directly match T1 cells to T0 cells using fuzzy lineage barcode overlap "
            "from Step 9 corrected UMI TSV files."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io_group = parser.add_argument_group("Input/output")
    io_group.add_argument(
        "--t0",
        required=True,
        help="Path to T0 Step 9 file (e.g., T0_corrected_UMI_counts.tsv)",
    )
    io_group.add_argument(
        "--t1-ctrl",
        default=None,
        help="Path to T1 control Step 9 TSV file",
    )
    io_group.add_argument(
        "--t1-treatment",
        default=None,
        help="Path to T1 treatment Step 9 TSV file",
    )
    io_group.add_argument(
        "--output-dir",
        required=True,
        help="Directory where output CSV files will be written",
    )

    matching_group = parser.add_argument_group("Matching parameters")
    matching_group.add_argument(
        "--hamming-threshold",
        type=non_negative_int,
        default=2,
        help="Maximum Hamming distance allowed for fuzzy barcode match",
    )
    matching_group.add_argument(
        "--jaccard-threshold",
        type=bounded_float_0_1,
        default=0.6,
        help="Minimum similarity required for a T0/T1 match to be retained",
    )
    matching_group.add_argument(
        "--min-matched-barcodes",
        type=positive_int,
        default=2,
        help="Minimum number of matched barcodes required for a match",
    )
    matching_group.add_argument(
        "--top-n",
        type=non_negative_int,
        default=0,
        help=(
            "Maximum matches to keep per T1 cell. Use 0 to keep all passing matches."
        ),
    )
    matching_group.add_argument(
        "--use-umi-weighting",
        action="store_true",
        help=(
            "Use UMI-weighted similarity as threshold/ranking score. "
            "Standard Jaccard is still reported."
        ),
    )

    perf_group = parser.add_argument_group("Performance/logging")
    perf_group.add_argument(
        "--n-jobs",
        type=non_negative_int,
        default=1,
        help="Number of worker processes. Use 0 for auto (cpu_count - 1).",
    )
    perf_group.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity level",
    )

    args = parser.parse_args(argv)

    if not args.t1_ctrl and not args.t1_treatment:
        parser.error("At least one of --t1-ctrl or --t1-treatment must be provided.")

    return args


def validate_runtime_inputs(args: argparse.Namespace) -> None:
    """Validate runtime configuration and input files before execution."""
    paths_to_check = [Path(args.t0)]
    if args.t1_ctrl:
        paths_to_check.append(Path(args.t1_ctrl))
    if args.t1_treatment:
        paths_to_check.append(Path(args.t1_treatment))

    for path in paths_to_check:
        if not path.exists():
            raise FileNotFoundError(f"Input file does not exist: {path}")
        if path.stat().st_size == 0:
            logging.warning("Input file appears empty: %s", path)

    if args.n_jobs == 0:
        args.n_jobs = max(1, (os.cpu_count() or 2) - 1)

    if args.use_umi_weighting and args.hamming_threshold < 0:
        raise ValueError("--hamming-threshold must be >= 0")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Script entrypoint."""
    args = parse_args(argv)
    setup_logging(args.log_level)

    try:
        validate_runtime_inputs(args)

        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logging.info("Loading T0 data: %s", args.t0)
        t0_data = parse_step9_tsv(Path(args.t0))
        logging.info("Loaded T0 cells: %d", len(t0_data.cell_to_barcodes))

        comparisons = []
        if args.t1_ctrl:
            comparisons.append(("T1_Ctrl", Path(args.t1_ctrl), output_dir / "T0_to_T1_Ctrl_direct_matches.csv"))
        if args.t1_treatment:
            comparisons.append(("T1_Treatment", Path(args.t1_treatment), output_dir / "T0_to_T1_Treatment_direct_matches.csv"))

        for label, t1_path, output_csv in comparisons:
            logging.info("Loading %s data: %s", label, t1_path)
            t1_data = parse_step9_tsv(t1_path)
            logging.info("Loaded %s cells: %d", label, len(t1_data.cell_to_barcodes))

            _run_condition_matching(
                condition_label=label,
                t0_data=t0_data,
                t1_data=t1_data,
                output_path=output_csv,
                hamming_threshold=args.hamming_threshold,
                jaccard_threshold=args.jaccard_threshold,
                min_matched_barcodes=args.min_matched_barcodes,
                use_umi_weighting=args.use_umi_weighting,
                n_jobs=args.n_jobs,
                top_n=args.top_n,
            )

        logging.info("All matching tasks completed successfully.")
        return 0

    except Exception as exc:
        logging.exception("Fatal error: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
