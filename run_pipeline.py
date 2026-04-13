#!/usr/bin/env python3
"""
RewindPipeline Master Orchestrator
===================================
Runs all 10 steps of the REWIND lineage tracing pipeline sequentially.

This script chains together every processing step—from raw FASTQ extraction
through clone ID assignment—with configurable parameters, automatic file-path
wiring, progress updates, and error handling.

Usage
-----
    python run_pipeline.py \
        --fastq1 sample_R1.fastq.gz \
        --fastq2 sample_R2.fastq.gz \
        --instrument-run-flowcell-ID "INSTRUMENT:RUN:FLOWCELL:" \
        --technology 10xv3 \
        --whitelist-file barcode_translation.tsv.gz \
        --edrops-file edrops_filtered_cells.tsv \
        --starcode-path /path/to/starcode \
        --output-dir results/ \
        --sample-name my_sample

For high-MOI experiments (multiple integrations per cell), set --top-n to
match the expected MOI (e.g., --top-n 6).
"""

import argparse
import os
import subprocess
import sys
import time
import textwrap


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(step_num, total_steps, msg):
    """Print a timestamped progress message."""
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [{step_num}/{total_steps}] {msg}", flush=True)


def run_cmd(cmd, step_num, total_steps, description):
    """Run a shell command with error handling."""
    log(step_num, total_steps, f"START  — {description}")
    log(step_num, total_steps, f"  CMD: {cmd}")
    t0 = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.stdout.strip():
        for line in result.stdout.strip().splitlines():
            print(f"        stdout> {line}")
    if result.stderr.strip():
        for line in result.stderr.strip().splitlines():
            print(f"        stderr> {line}")

    if result.returncode != 0:
        print(f"\n{'='*60}")
        print(f"ERROR in step {step_num}: {description}")
        print(f"Return code: {result.returncode}")
        print(f"{'='*60}\n")
        sys.exit(result.returncode)

    log(step_num, total_steps, f"DONE   — {description} ({elapsed:.1f}s)")
    return result


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            RewindPipeline Master Orchestrator
            ==================================
            Runs all 10 pipeline steps sequentially with automatic
            file-path wiring and progress reporting.
        """),
    )

    # -- Required inputs --
    req = parser.add_argument_group("Required inputs")
    req.add_argument("--fastq1", required=True, help="Path to Read 1 FASTQ (.fastq.gz)")
    req.add_argument("--fastq2", required=True, help="Path to Read 2 FASTQ (.fastq.gz)")
    req.add_argument("--instrument-run-flowcell-ID", required=True,
                     help="Instrument:Run:Flowcell ID prefix from FASTQ headers")
    req.add_argument("--technology", required=True,
                     choices=["DropSeqv1", "DropSeqv2", "10xv2", "10xv3",
                              "CiteSeq5v2", "CiteSeq3v3", "CiteSeqTSB", "PearSeq"],
                     help="Library chemistry")
    req.add_argument("--whitelist-file", required=True,
                     help="Barcode translation whitelist (TSV, optionally .gz)")
    req.add_argument("--edrops-file", required=True,
                     help="eDROPs-filtered cell barcode list (TSV)")
    req.add_argument("--starcode-path", required=True,
                     help="Path to compiled Starcode executable")

    # -- Output --
    out = parser.add_argument_group("Output")
    out.add_argument("--output-dir", default="rewind_output",
                     help="Directory for all intermediate and final outputs (default: rewind_output)")
    out.add_argument("--sample-name", default="rewind",
                     help="Prefix for output file names (default: rewind)")

    # -- Tunable parameters --
    tune = parser.add_argument_group("Tunable parameters")
    tune.add_argument("--capture-seq", default="GCTCACCTATTAGCGGCTAAGG",
                      help="Capture sequence used to locate lineage barcodes in Read 2 "
                           "(default: GCTCACCTATTAGCGGCTAAGG)")
    tune.add_argument("--barcode-length", type=int, default=84,
                      help="Length of the lineage barcode (nt) upstream of the capture "
                           "sequence (default: 84)")
    tune.add_argument("--max-mismatches", type=int, default=0,
                      help="Max WSN mismatches for Step 2 filtering (default: 0)")
    tune.add_argument("--min-reads-count", type=int, default=1,
                      help="Min reads per molecule for Step 6 (default: 1)")
    tune.add_argument("--starcode-max-distance", type=int, default=3,
                      help="Levenshtein distance for Starcode clustering (default: 3)")
    tune.add_argument("--top-n", type=int, default=1,
                      help="Number of top barcodes to estimate expected count in "
                           "noise filtering. Set to ~MOI for high-MOI experiments "
                           "(default: 1)")
    tune.add_argument("--z-threshold", type=float, default=2.0,
                      help="Z-score threshold for noise filtering (default: 2.0)")
    tune.add_argument("--min-group-size", type=int, default=5,
                      help="Minimum subpopulation size for Z-score calculation (default: 5)")
    tune.add_argument("--max-lineage-count", type=int, default=None,
                      help="Max lineage barcodes per cell for clone assignment (optional)")

    # -- Timepoint matching (optional Step 11) --
    tp = parser.add_argument_group("Timepoint matching (optional)")
    tp.add_argument("--run-timepoint-matching", action="store_true",
                    help="Run timepoint matching after Step 10 (requires --t0-reference)")
    tp.add_argument("--t0-reference", default=None,
                    help="Path to T0 cloneID_cloneBarcode.tsv for timepoint matching")
    tp.add_argument("--t0-filtered-reference", default=None,
                    help="(Optional) T0 rewind_filtered.tsv for UMI-weighted matching")
    tp.add_argument("--hamming-threshold", type=int, default=2,
                    help="Hamming distance threshold for fuzzy barcode matching (default: 2)")
    tp.add_argument("--similarity-threshold", type=float, default=0.8,
                    help="Minimum similarity for a valid match (default: 0.8)")
    tp.add_argument("--top-n-matches", type=int, default=5,
                    help="Top-N ranked matches per T1 clone (default: 5)")
    tp.add_argument("--weight-by-umi", action="store_true",
                    help="Weight similarity by UMI counts")
    tp.add_argument("--matching-n-jobs", type=int, default=1,
                    help="Parallel workers for timepoint matching (default: 1; 0=auto)")

    # -- Differential enrichment analysis (optional Step 12) --
    ea = parser.add_argument_group("Differential enrichment analysis (optional)")
    ea.add_argument("--run-enrichment-analysis", action="store_true",
                    help="Run differential enrichment analysis after timepoint matching "
                         "(requires --control-reference and --treated-reference)")
    ea.add_argument("--control-reference", default=None,
                    help="Path to T1 Control cloneID_cloneBarcode.tsv file")
    ea.add_argument("--treated-reference", default=None,
                    help="Path to T1 Treated cloneID_cloneBarcode.tsv file "
                         "(omit if the current pipeline run IS the treated sample)")
    ea.add_argument("--control-matches", default=None,
                    help="Path to Control best_matches.txt from timepoint matching")
    ea.add_argument("--treated-matches-file", default=None,
                    help="Path to Treated best_matches.txt from timepoint matching")
    ea.add_argument("--control-filtered-reference", default=None,
                    help="Path to Control rewind_filtered.tsv for UMI-based sizing")
    ea.add_argument("--treated-filtered-reference", default=None,
                    help="Path to Treated rewind_filtered.tsv for UMI-based sizing")
    ea.add_argument("--pseudocount", type=float, default=1.0,
                    help="Pseudocount for fold-change calculation (default: 1)")
    ea.add_argument("--fdr-threshold", type=float, default=0.05,
                    help="FDR q-value threshold for significance (default: 0.05)")
    ea.add_argument("--log2fc-threshold", type=float, default=1.0,
                    help="Log2 fold-change ratio threshold (default: 1.0)")

    # -- Control --
    ctrl = parser.add_argument_group("Pipeline control")
    ctrl.add_argument("--start-step", type=int, default=1, choices=range(1, 13),
                      help="Step to start from (1–12). Useful for resuming. "
                           "Intermediate files must already exist. (default: 1)")
    ctrl.add_argument("--stop-step", type=int, default=10, choices=range(1, 13),
                      help="Step to stop after (1–12; 11=timepoint matching, "
                           "12=enrichment analysis). (default: 10)")
    ctrl.add_argument("--python", default=sys.executable,
                      help="Python interpreter to use (default: current interpreter)")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    TOTAL = 10
    py = args.python
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Derive intermediate file names
    s = args.sample_name
    d = args.output_dir

    f_extracted       = os.path.join(d, f"{s}_address.tsv")
    f_filtered        = os.path.join(d, f"{s}_address.mm{args.max_mismatches}.tsv")
    f_fcts            = os.path.join(d, f"{s}_address.mm{args.max_mismatches}.fcts.txt")
    f_collapsed       = os.path.join(d, f"{s}_address.mm{args.max_mismatches}.fcts.collapsed.txt")
    f_trackbcs        = os.path.join(d, f"{s}_address.mm{args.max_mismatches}.fcts.trackbcs.txt")
    f_recollapsed     = os.path.join(d, f"{s}_address.mm{args.max_mismatches}.fcts.collapsed.recollapsed.txt")
    f_translated      = os.path.join(d, f"{s}_translated_filtered.tsv")
    f_umi_counts      = os.path.join(d, f"{s}_UMI_counts.txt")
    f_starcode_out    = os.path.join(d, f"{s}_starcode.txt")
    f_corrected       = os.path.join(d, f"{s}_corrected_UMI_counts.tsv")
    # Step 10 outputs are derived from f_corrected by the script itself

    print(f"\n{'='*60}")
    print(f"  RewindPipeline — Master Orchestrator")
    print(f"  Sample:     {s}")
    print(f"  Output dir: {os.path.abspath(d)}")
    print(f"  Steps:      {args.start_step} → {args.stop_step}")
    print(f"  top-n:      {args.top_n}" + (" (HIGH MOI MODE)" if args.top_n > 1 else ""))
    print(f"{'='*60}\n")

    t_start = time.time()

    # --- Step 1 ---
    if args.start_step <= 1 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "rewind_extractor.py")} '
            f'--fastq1 {args.fastq1} '
            f'--fastq2 {args.fastq2} '
            f'--output {f_extracted} '
            f'--instrument-run-flowcell-ID "{args.instrument_run_flowcell_ID}" '
            f'--capture-seq {args.capture_seq} '
            f'--barcode-length {args.barcode_length}',
            1, TOTAL, "Extract barcodes from FASTQ files"
        )

    # --- Step 2 ---
    if args.start_step <= 2 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "filter_by_mismatch.py")} '
            f'--input {f_extracted} '
            f'--output {f_filtered} '
            f'--max_mismatches {args.max_mismatches}',
            2, TOTAL, "Filter by WSN pattern mismatches"
        )

    # --- Step 3 ---
    if args.start_step <= 3 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "addressct_rewind.py")} '
            f'--input {f_filtered} '
            f'--output {f_fcts}',
            3, TOTAL, "First-pass UMI collapse"
        )

    # --- Step 4 ---
    if args.start_step <= 4 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "collapse.py")} '
            f'{f_fcts} {f_collapsed} {f_trackbcs} {args.technology}',
            4, TOTAL, "Correct cell barcodes (incomplete extension)"
        )

    # --- Step 5 ---
    if args.start_step <= 5 <= args.stop_step:
        run_cmd(
            f'cat {f_collapsed} | {py} {os.path.join(pipeline_dir, "addressct2.py")} {f_recollapsed}',
            5, TOTAL, "Second-pass UMI collapse"
        )

    # --- Step 6 ---
    if args.start_step <= 6 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "translate_and_filter_cell_barcodes.py")} '
            f'--rewind_file {f_recollapsed} '
            f'--whitelist_file {args.whitelist_file} '
            f'--edrops_file {args.edrops_file} '
            f'--output_file {f_translated} '
            f'--min-reads-count {args.min_reads_count}',
            6, TOTAL, "Translate & filter cell barcodes"
        )

    # --- Step 7 ---
    if args.start_step <= 7 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "rewind_address_UMI_counts.py")} '
            f'--input {f_translated} '
            f'--output {f_umi_counts}',
            7, TOTAL, "Summarize UMI counts per cell–lineage pair"
        )

    # --- Step 8 ---
    if args.start_step <= 8 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "run_starcode_on_rewind.py")} '
            f'--input {f_umi_counts} '
            f'--output {f_starcode_out} '
            f'--starcode_path {args.starcode_path} '
            f'--max_distance {args.starcode_max_distance}',
            8, TOTAL, "Error-correct lineage barcodes with Starcode"
        )

    # --- Step 9 ---
    if args.start_step <= 9 <= args.stop_step:
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "correct_rewind_and_collapse_umi.py")} '
            f'--starcode_file {f_starcode_out} '
            f'--umi_counts_file {f_umi_counts} '
            f'--output_file {f_corrected}',
            9, TOTAL, "Apply Starcode corrections & collapse UMIs"
        )

    # --- Step 10 ---
    if args.start_step <= 10 <= args.stop_step:
        mlc_flag = f"--max-lineage-count {args.max_lineage_count}" if args.max_lineage_count else ""
        run_cmd(
            f'{py} {os.path.join(pipeline_dir, "rewind_barcode_filtering_and_clone_assignment.py")} '
            f'--input-file {f_corrected} '
            f'--top-n {args.top_n} '
            f'--z-threshold {args.z_threshold} '
            f'--min-group-size {args.min_group_size} '
            f'{mlc_flag}',
            10, TOTAL, "Filter noise & assign clone IDs"
        )

    # --- Step 11 (optional): Timepoint matching ---
    if args.run_timepoint_matching and args.start_step <= 11 <= args.stop_step:
        if not args.t0_reference:
            print("WARNING: --run-timepoint-matching requires --t0-reference. Skipping Step 11.")
        else:
            # Derive T1 clone barcode file from Step 10 outputs
            corrected_base = os.path.splitext(f_corrected)[0]
            z_suffix = f"z{args.z_threshold}"
            mlc_suffix = f"_mlc{args.max_lineage_count}" if args.max_lineage_count else ""
            t1_clone_barcode = f"{corrected_base}.{z_suffix}{mlc_suffix}.cloneID_cloneBarcode.tsv"
            t1_filtered = f"{corrected_base}.{z_suffix}.rewind_filtered.tsv"

            matching_output = os.path.join(d, "timepoint_matching")

            umi_flags = ""
            if args.weight_by_umi:
                t0_filt = args.t0_filtered_reference or ""
                if t0_filt:
                    umi_flags = (
                        f"--weight-by-umi "
                        f"--t0-filtered {t0_filt} "
                        f"--t1-filtered {t1_filtered} "
                    )

            run_cmd(
                f'{py} {os.path.join(pipeline_dir, "match_timepoints.py")} '
                f'--t0-clone-barcode {args.t0_reference} '
                f'--t1-clone-barcode {t1_clone_barcode} '
                f'--output-dir {matching_output} '
                f'--hamming-threshold {args.hamming_threshold} '
                f'--similarity-threshold {args.similarity_threshold} '
                f'--top-n-matches {args.top_n_matches} '
                f'--n-jobs {args.matching_n_jobs} '
                f'{umi_flags}',
                11, TOTAL + 1, "Timepoint matching (T1 → T0)"
            )

    # --- Step 12 (optional): Differential enrichment analysis ---
    if args.run_enrichment_analysis and args.start_step <= 12 <= args.stop_step:
        if not args.t0_reference:
            print("WARNING: --run-enrichment-analysis requires --t0-reference. Skipping Step 12.")
        elif not args.control_reference:
            print("WARNING: --run-enrichment-analysis requires --control-reference. Skipping Step 12.")
        else:
            # Derive file paths
            corrected_base = os.path.splitext(f_corrected)[0]
            z_suffix = f"z{args.z_threshold}"
            mlc_suffix = f"_mlc{args.max_lineage_count}" if args.max_lineage_count else ""

            # T0 clone barcode file is the user-supplied --t0-reference
            t0_clones = args.t0_reference

            # Control clone barcode file
            ctrl_clones = args.control_reference

            # Treated clone barcode file: use --treated-reference, or this run's own output
            if args.treated_reference:
                treat_clones = args.treated_reference
            else:
                treat_clones = f"{corrected_base}.{z_suffix}{mlc_suffix}.cloneID_cloneBarcode.tsv"

            # Matching results
            matching_output = os.path.join(d, "timepoint_matching")
            ctrl_matches = args.control_matches or os.path.join(matching_output, "control_best_matches.txt")
            treat_matches = args.treated_matches_file or os.path.join(matching_output, "best_matches.txt")

            enrichment_output = os.path.join(d, "enrichment_analysis")

            # Build optional filtered-file flags
            filt_flags = ""
            t0_filt = args.t0_filtered_reference
            ctrl_filt = args.control_filtered_reference
            treat_filt = args.treated_filtered_reference
            if not treat_filt:
                # If current run is the treated sample, use its own filtered file
                treat_filt_auto = f"{corrected_base}.{z_suffix}.rewind_filtered.tsv"
                if os.path.isfile(treat_filt_auto):
                    treat_filt = treat_filt_auto
            if t0_filt:
                filt_flags += f"--t0-filtered {t0_filt} "
            if ctrl_filt:
                filt_flags += f"--control-filtered {ctrl_filt} "
            if treat_filt:
                filt_flags += f"--treated-filtered {treat_filt} "

            run_cmd(
                f'{py} {os.path.join(pipeline_dir, "differential_enrichment_analysis.py")} '
                f'--t0-clones {t0_clones} '
                f'--control-clones {ctrl_clones} '
                f'--treated-clones {treat_clones} '
                f'--control-matches {ctrl_matches} '
                f'--treated-matches {treat_matches} '
                f'--output-dir {enrichment_output} '
                f'--pseudocount {args.pseudocount} '
                f'--fdr-threshold {args.fdr_threshold} '
                f'--log2fc-threshold {args.log2fc_threshold} '
                f'{filt_flags}',
                12, TOTAL + 2, "Differential enrichment analysis"
            )

    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"  Pipeline complete!  Total time: {elapsed:.1f}s")
    print(f"  Outputs in: {os.path.abspath(d)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
