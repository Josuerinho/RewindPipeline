# RewindPipeline

A computational pipeline for processing **lineage tracing single-cell RNA-seq (scRNA-seq) data** generated using the **REWIND** (RNA Expression and DNA barcode recovery With IN situ Dissection) method. The pipeline extracts lineage barcodes from FASTQ files, performs quality filtering, error correction, UMI collapsing, and ultimately assigns cells to clonal identities.

---

## Overview

REWIND is a lineage tracing technique that uses lentiviral-delivered DNA barcodes to track clonal relationships between single cells. In a typical experiment, cells are transduced with a barcoded lentiviral library, expanded, and then profiled by scRNA-seq (e.g., 10x Genomics or Drop-seq). This pipeline processes the resulting sequencing data to:

1. **Extract** lineage barcodes and cell barcodes from raw FASTQ files
2. **Filter** reads by barcode quality (WSN pattern mismatches)
3. **Collapse** UMIs to count unique molecules per cell–barcode pair
4. **Correct** cell barcodes for synthesis errors (incomplete extension)
5. **Translate** feature barcoding cell barcodes to gene expression cell barcodes
6. **Error-correct** lineage barcodes using Starcode clustering
7. **Filter noise** and **assign clone IDs** to cells

---

## Pipeline Execution Order

The pipeline consists of scripts that should be run **sequentially** in the following order. You can run them individually via the command line, or use the **`run_pipeline.py` master orchestrator** (see below) to run all steps automatically.

```
Step 1: rewind_extractor.py           → Extract barcodes from FASTQ files
Step 2: filter_by_mismatch.py         → Filter by WSN pattern mismatches
Step 3: addressct_rewind.py           → Collapse read counts by UMI (first pass)
Step 4: collapse.py                   → Correct cell barcodes (incomplete extension)
Step 5: addressct2.py                 → Re-collapse UMIs after barcode correction
Step 6: translate_and_filter_cell_barcodes.py → Translate & filter cell barcodes
Step 7: rewind_address_UMI_counts.py  → Summarize UMI counts per cell–lineage pair
Step 8: run_starcode_on_rewind.py     → Error-correct lineage barcodes with Starcode
Step 9: correct_rewind_and_collapse_umi.py → Apply Starcode corrections & collapse
Step 10: rewind_barcode_filtering_and_clone_assignment.py → Filter noise & assign clone IDs
```

Utility/QC scripts (run as needed):
- `compute_reads_percentage.py` — Compute percentage of reads retained after a processing step
- `dsstats.py` — Generate summary statistics and plots for a Drop-seq experiment

---

## Repository Structure

```
RewindPipeline/
├── rewind_extractor.py                          # Step 1: Extract barcodes from FASTQ
├── filter_by_mismatch.py                        # Step 2: Filter by mismatch count
├── addressct_rewind.py                          # Step 3: First UMI collapse
├── collapse.py                                  # Step 4: Cell barcode correction
├── addressct2.py                                # Step 5: Second UMI collapse
├── translate_and_filter_cell_barcodes.py        # Step 6: Barcode translation & filtering
├── rewind_address_UMI_counts.py                 # Step 7: UMI count summary
├── run_starcode_on_rewind.py                    # Step 8: Starcode lineage barcode correction
├── correct_rewind_and_collapse_umi.py           # Step 9: Apply corrections & collapse UMIs
├── rewind_barcode_filtering_and_clone_assignment.py  # Step 10: Noise filtering & clone assignment
├── run_pipeline.py                                   # Master orchestrator (runs all 10 steps)
├── compute_reads_percentage.py                  # Utility: Read retention statistics
├── dsstats.py                                   # Utility: Drop-seq experiment QC plots
└── .gitignore
```

---

## Detailed Script Descriptions

### Step 1: `rewind_extractor.py`
**Purpose:** Extracts cell barcodes (16 nt), UMIs (12 nt), and 84 nt lineage barcodes from paired-end FASTQ files (10x Genomics format). Locates lineage barcodes in Read 2 by searching for the capture sequence (`GCTCACCTATTAGCGGCTAAGG`) and extracting the 84 nt upstream. Counts mismatches to the expected **WSN trinucleotide repeat pattern** (W = A/T, S = G/C, N = any).

```bash
python rewind_extractor.py \
  --fastq1 <read1.fastq.gz> \
  --fastq2 <read2.fastq.gz> \
  --output <rewind_address.tsv> \
  --instrument-run-flowcell-ID <e.g., AV100007:PY055:2336402118:>
```

**Output format:** `readid \t cell_barcode \t umi \t lineage_barcode \t mismatches`

### Step 2: `filter_by_mismatch.py`
**Purpose:** Filters extracted reads, keeping only those where the lineage barcode has ≤ N mismatches to the WSN pattern.

```bash
python filter_by_mismatch.py \
  --input <rewind_address.tsv> \
  --output <rewind_address.filtered.tsv> \
  --max_mismatches <0>
```

### Step 3: `addressct_rewind.py`
**Purpose:** First-pass UMI collapsing — counts the number of reads per unique (cell_barcode, UMI, lineage_barcode) combination. Supports `.gz`, `.bz2`, and plain text input.

```bash
python addressct_rewind.py \
  --input <rewind_address.filtered.tsv> \
  --output <rewind_address.fcts.txt>
```

**Output format:** `cell_barcode \t umi \t lineage_barcode \t read_count`

### Step 4: `collapse.py`
**Purpose:** Identifies and corrects cell barcodes affected by **incomplete oligonucleotide extension** during synthesis. Detects barcodes where >90% of UMI terminal nucleotides are 'T' (indicating the cell barcode is 1 nt too short, with the missing nt replaced by oligo-dT priming). Collapses the four possible erroneous 12/16 nt barcodes into the correct shorter barcode and adjusts UMI sequences accordingly.

Supports multiple library chemistries: `DropSeqv1`, `DropSeqv2`, `10xv2`, `10xv3`, `CiteSeq5v2`, `CiteSeq3v3`, `CiteSeqTSB`, `PearSeq`.

```bash
python collapse.py \
  <input_addresscts.txt> \
  <output_collapsed.txt> \
  <output_trackbcs.txt> \
  <technology>  # e.g., 10xv3
```

### Step 5: `addressct2.py`
**Purpose:** Second-pass UMI collapsing after cell barcode correction. Since `collapse.py` modifies cell barcode and UMI sequences, this re-collapses identical (cell_barcode, UMI, gene/lineage) addresses. Reads from stdin.

```bash
cat <collapsed_output.txt> | python addressct2.py <output_recollapsed.txt>
```

### Step 6: `translate_and_filter_cell_barcodes.py`
**Purpose:** Translates cell barcodes from the **feature barcoding** library to the corresponding **gene expression** library barcodes using a whitelist mapping file. Then filters to retain only cell barcodes present in an eDROPs-filtered cell list (i.e., cells that pass quality filtering in the gene expression data).

```bash
python translate_and_filter_cell_barcodes.py \
  --rewind_file <recollapsed_address.txt> \
  --whitelist_file <barcode_translation_whitelist.tsv[.gz]> \
  --edrops_file <edrops_filtered_cells.tsv> \
  --output_file <translated_filtered.tsv> \
  --min-reads-count <1>
```

### Step 7: `rewind_address_UMI_counts.py`
**Purpose:** Summarizes the data by counting the number of **unique UMIs** for each (cell_barcode, lineage_barcode) pair. Also computes total unique UMIs and total unique lineage barcodes per cell. Output is sorted by cells with the most molecules.

```bash
python rewind_address_UMI_counts.py \
  --input <translated_filtered.tsv> \
  --output <UMI_counts.txt>
```

### Step 8: `run_starcode_on_rewind.py`
**Purpose:** Error-corrects lineage barcodes within each cell barcode using **[Starcode](https://github.com/gui11aume/starcode)**, a sequence clustering tool based on Levenshtein distance. Groups similar lineage barcode sequences and maps each to a canonical (consensus) sequence.

```bash
python run_starcode_on_rewind.py \
  --input <UMI_counts.txt> \
  --output <starcode_output.txt> \
  --starcode_path </path/to/starcode> \
  --max_distance <3>
```

**Requires:** The [Starcode](https://github.com/gui11aume/starcode) executable installed and accessible.

### Step 9: `correct_rewind_and_collapse_umi.py`
**Purpose:** Applies the Starcode error corrections to lineage barcodes in the UMI counts file. Merges UMI counts for barcodes that were clustered together, then computes summary statistics (total molecules and lineage barcodes per cell).

```bash
python correct_rewind_and_collapse_umi.py \
  --starcode_file <starcode_output.txt> \
  --umi_counts_file <UMI_counts.txt> \
  --output_file <corrected_UMI_counts.tsv>
```

### Step 10: `rewind_barcode_filtering_and_clone_assignment.py`
**Purpose:** Final analysis step that performs two major tasks:

1. **Noise filtering:** Normalizes UMI counts per cell, log-transforms, and calculates Z-scores to identify and remove low-confidence lineage barcode assignments. Tests if the resulting distribution of lineage barcodes per cell fits a **truncated Poisson distribution** (expected for random lentiviral integration).

2. **Clone ID assignment:** Determines clone barcodes for each cell (concatenating sorted lineage barcodes for cells with multiple integrations), and assigns a unique clone ID to each distinct clone barcode. Generates QC plots.

```bash
python rewind_barcode_filtering_and_clone_assignment.py \
  --input-file <corrected_UMI_counts.tsv> \
  --top-n 1 \
  --z-threshold 2 \
  --min-group-size 5 \
  --max-lineage-count <optional_max>
```

**Outputs** (using input file basename as prefix):
- `*.rewind_statistics.tsv` — Full statistics for all lineage barcode observations
- `*.z{threshold}.rewind_filtered.tsv` — Filtered lineage barcode data
- `*.z{threshold}.cloneID_cloneBarcode.tsv` — Clone ID to clone barcode mapping
- `*.z{threshold}.cloneID_cellBarcode.tsv` — Clone ID to cell barcode mapping
- `*.truncated_Poisson_fit_results.txt` — Poisson fit QC statistics
- `*.truncated_Poisson_fit_plots.png` — Poisson fit QC plots
- `*.clone_barcode_distribution_qc_plot.png` — Clone size distribution plot

### Utility: `compute_reads_percentage.py`
**Purpose:** Computes the percentage of reads retained after a processing step, comparing an original file (FASTQ or text) to a processed file.

```bash
python compute_reads_percentage.py \
  --original <original_file.fastq.gz> \
  --processed <processed_file.tsv>
```

### Utility: `dsstats.py`
**Purpose:** Generates summary statistics and PDF plots (cumulative molecule histograms, molecules/genes per cell barcode) for a Drop-seq experiment.

```bash
python dsstats.py <filtered_addresses.txt> <output.pdf> <cumulative_hist.txt>
```

---

## Dependencies

### Python Packages

| Package | Used By |
|---------|---------|
| `numpy` | `collapse.py`, `rewind_barcode_filtering_and_clone_assignment.py`, `correct_rewind_and_collapse_umi.py` |
| `pandas` | `translate_and_filter_cell_barcodes.py`, `rewind_address_UMI_counts.py`, `run_starcode_on_rewind.py`, `correct_rewind_and_collapse_umi.py`, `rewind_barcode_filtering_and_clone_assignment.py` |
| `biopython` (`Bio.SeqIO`) | `rewind_extractor.py` |
| `matplotlib` | `dsstats.py`, `rewind_barcode_filtering_and_clone_assignment.py` |
| `seaborn` | `rewind_barcode_filtering_and_clone_assignment.py` |
| `scipy` | `rewind_barcode_filtering_and_clone_assignment.py` |

Install all dependencies:

```bash
pip install numpy pandas biopython matplotlib seaborn scipy
```

### External Tools

- **[Starcode](https://github.com/gui11aume/starcode)** — Required for Step 8 (`run_starcode_on_rewind.py`). A sequence clustering tool that groups similar sequences using Levenshtein distance. Must be compiled and accessible as an executable.

  ```bash
  git clone https://github.com/gui11aume/starcode.git
  cd starcode
  make
  ```

---

## Input Data Requirements

| Input | Description |
|-------|-------------|
| **Read 1 FASTQ** (`.fastq.gz`) | Contains cell barcodes (first 16 nt) and UMIs (next 12 nt) |
| **Read 2 FASTQ** (`.fastq.gz`) | Contains lineage barcodes (84 nt upstream of capture sequence) |
| **Instrument/Run/Flowcell ID** | Header prefix from FASTQ read IDs (e.g., `AV100007:PY055:2336402118:`) |
| **Barcode translation whitelist** | TSV mapping feature barcoding cell barcodes to gene expression cell barcodes (plain or `.gz`) |
| **eDROPs filtered cell list** | TSV of cell barcodes that pass gene expression QC (barcode, molecule count, gene count) |
| **Starcode executable** | Path to compiled Starcode binary |

---

## Quick Start Example

```bash
# Step 1: Extract barcodes from FASTQ files
python rewind_extractor.py \
  --fastq1 sample_R1.fastq.gz \
  --fastq2 sample_R2.fastq.gz \
  --output rewind_address.tsv \
  --instrument-run-flowcell-ID "INSTRUMENT:RUN:FLOWCELL:"

# Step 2: Filter by WSN mismatch (keep perfect matches only)
python filter_by_mismatch.py \
  --input rewind_address.tsv \
  --output rewind_address.mm0.tsv \
  --max_mismatches 0

# Step 3: First UMI collapse
python addressct_rewind.py \
  --input rewind_address.mm0.tsv \
  --output rewind_address.mm0.fcts.txt

# Step 4: Correct cell barcodes for incomplete extension
python collapse.py \
  rewind_address.mm0.fcts.txt \
  rewind_address.mm0.fcts.collapsed.txt \
  rewind_address.mm0.fcts.trackbcs.txt \
  10xv3

# Step 5: Second UMI collapse (after barcode correction)
cat rewind_address.mm0.fcts.collapsed.txt | python addressct2.py rewind_address.mm0.fcts.collapsed.recollapsed.txt

# Step 6: Translate feature barcodes → GEX barcodes, and filter
python translate_and_filter_cell_barcodes.py \
  --rewind_file rewind_address.mm0.fcts.collapsed.recollapsed.txt \
  --whitelist_file barcode_translation.tsv.gz \
  --edrops_file edrops_filtered_cells.tsv \
  --output_file rewind_translated_filtered.tsv \
  --min-reads-count 1

# Step 7: Summarize UMI counts per cell-lineage pair
python rewind_address_UMI_counts.py \
  --input rewind_translated_filtered.tsv \
  --output rewind_UMI_counts.txt

# Step 8: Error-correct lineage barcodes with Starcode
python run_starcode_on_rewind.py \
  --input rewind_UMI_counts.txt \
  --output rewind_starcode.txt \
  --starcode_path /path/to/starcode \
  --max_distance 3

# Step 9: Apply Starcode corrections
python correct_rewind_and_collapse_umi.py \
  --starcode_file rewind_starcode.txt \
  --umi_counts_file rewind_UMI_counts.txt \
  --output_file rewind_corrected.tsv

# Step 10: Filter noise and assign clone IDs
python rewind_barcode_filtering_and_clone_assignment.py \
  --input-file rewind_corrected.tsv \
  --top-n 1 \
  --z-threshold 2 \
  --min-group-size 5
```

---

## Key Parameters to Tune

| Parameter | Script | Default | Description |
|-----------|--------|---------|-------------|
| `--max_mismatches` | `filter_by_mismatch.py` | — | Max allowed WSN pattern mismatches (0 = perfect match only) |
| `technology` | `collapse.py` | — | Library chemistry: `10xv2`, `10xv3`, `DropSeqv1`, `DropSeqv2`, `PearSeq`, `CiteSeq5v2`, `CiteSeq3v3`, `CiteSeqTSB` |
| `--min-reads-count` | `translate_and_filter_cell_barcodes.py` | — | Minimum reads per molecule to retain |
| `--max_distance` | `run_starcode_on_rewind.py` | — | Levenshtein distance for Starcode clustering |
| `--z-threshold` | `rewind_barcode_filtering_and_clone_assignment.py` | 2 | Z-score threshold for noise filtering |
| `--top-n` | `rewind_barcode_filtering_and_clone_assignment.py` | 1 | Number of top barcodes used to estimate expected count |
| `--max-lineage-count` | `rewind_barcode_filtering_and_clone_assignment.py` | None | Max lineage barcodes per cell for clone assignment |

---

## Master Orchestrator Script (`run_pipeline.py`)

A convenience script that chains all 10 steps into a single command with automatic file-path wiring, progress reporting, and error handling.

### Basic usage

```bash
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
```

### High-MOI usage

```bash
python run_pipeline.py \
  --fastq1 sample_R1.fastq.gz \
  --fastq2 sample_R2.fastq.gz \
  --instrument-run-flowcell-ID "INSTRUMENT:RUN:FLOWCELL:" \
  --technology 10xv3 \
  --whitelist-file barcode_translation.tsv.gz \
  --edrops-file edrops_filtered_cells.tsv \
  --starcode-path /path/to/starcode \
  --output-dir results/ \
  --sample-name my_sample \
  --top-n 6 \
  --z-threshold 2 \
  --max-lineage-count 15
```

### Features

- **All parameters configurable** — every tunable parameter from the individual scripts is exposed as a CLI flag
- **Automatic file naming** — intermediate files are placed in `--output-dir` with `--sample-name` as prefix
- **Resumable** — use `--start-step N` and `--stop-step M` to run a subset of steps (intermediate files must exist)
- **Progress reporting** — timestamped log messages for each step with elapsed time
- **Error handling** — the pipeline stops immediately if any step fails, reporting the failing command

### Full parameter list

| Flag | Default | Description |
|------|---------|-------------|
| `--fastq1` | *required* | Path to Read 1 FASTQ |
| `--fastq2` | *required* | Path to Read 2 FASTQ |
| `--instrument-run-flowcell-ID` | *required* | FASTQ header prefix |
| `--technology` | *required* | Library chemistry |
| `--whitelist-file` | *required* | Barcode translation whitelist |
| `--edrops-file` | *required* | eDROPs-filtered cell list |
| `--starcode-path` | *required* | Path to Starcode binary |
| `--output-dir` | `rewind_output` | Output directory |
| `--sample-name` | `rewind` | File name prefix |
| `--max-mismatches` | `0` | WSN mismatch threshold |
| `--min-reads-count` | `1` | Min reads per molecule |
| `--starcode-max-distance` | `3` | Starcode Levenshtein distance |
| `--top-n` | `1` | Top barcodes for expected count (set to ~MOI for high-MOI) |
| `--z-threshold` | `2.0` | Z-score noise filter threshold |
| `--min-group-size` | `5` | Min subpopulation size for Z-score |
| `--max-lineage-count` | *None* | Max lineage count for clone assignment |
| `--start-step` | `1` | Step to start from (1–10) |
| `--stop-step` | `10` | Step to stop after (1–10) |
| `--python` | *current* | Python interpreter path |

### Automation note

All 10 steps are **fully automatable** — none require human-in-the-loop inspection or manual decision-making to proceed. However, it is **strongly recommended** to inspect intermediate QC outputs after the first run:

- After **Step 2**: Check read retention rates (use `compute_reads_percentage.py`)
- After **Step 9**: Review the fraction of cells with single vs. multiple lineage barcodes
- After **Step 10**: Inspect the truncated Poisson fit plots and clone size distribution

These QC checks help validate parameter choices but do not block pipeline execution.

---

## High MOI (Multiple Integrations per Cell) Scenarios

The default pipeline configuration assumes **single lentiviral integration** per cell (MOI ≈ 1), using `--top-n 1` in the final noise-filtering step. For experiments with **high MOI** (e.g., 6–7 integrations per cell), the pipeline can be adapted with the following guidance.

### What `--top-n` controls

In Step 10 (`rewind_barcode_filtering_and_clone_assignment.py`), the noise filter works as follows:

1. **Normalize** UMI counts per cell to a fixed target sum (10,000)
2. **Log-transform** the normalized counts
3. **Estimate expected count** per cell: the mean of the **top N** log-normalized lineage barcode counts for that cell
4. **Calculate Z-scores** to identify lineage barcodes that deviate significantly from the expected count
5. **Filter** entries where |Z-score| ≥ threshold

With `--top-n 1`, only the single highest-count barcode is used to estimate what a "real" barcode looks like. For high-MOI experiments, you should set `--top-n` to approximately match your expected MOI so that multiple real barcodes contribute to the expected count estimate.

### Recommended parameter changes for high MOI

| Parameter | Low MOI (default) | High MOI (~6–7) | Rationale |
|-----------|-------------------|------------------|-----------|
| `--top-n` | 1 | 6 or 7 | Match your expected MOI so all real barcodes contribute to the expected count |
| `--z-threshold` | 2.0 | 2.0–3.0 | You may want a slightly more permissive threshold since the variance across real barcodes is higher |
| `--max-lineage-count` | None | 15–20 | Exclude cells with implausibly many lineage barcodes (e.g., doublets or noise); set to ~2–3× your expected MOI |

### How the pipeline handles multiple barcodes per cell

The pipeline is architecturally compatible with multiple barcodes per cell throughout all 10 steps:

- **Steps 1–9** (extraction through error correction) are barcode-count agnostic — they process every (cell, lineage_barcode) pair independently and do not assume one barcode per cell.
- **Step 10 — Noise filtering**: With `--top-n` set appropriately, the Z-score calculation correctly identifies real vs. noise barcodes even when cells have multiple real integrations.
- **Step 10 — Clone ID assignment**: The `generate_clone_barcodes()` function already handles multiple barcodes per cell by concatenating sorted lineage barcodes with `_` as separator (e.g., `BARCODE_A_BARCODE_B_BARCODE_C`). Cells sharing the same set of barcodes are assigned the same clone ID.

### Truncated Poisson model at high MOI

The truncated Poisson fit (QC check in Step 10) tests whether the observed distribution of lineage barcodes per cell fits a zero-truncated Poisson distribution. This model:

- **Still applies at high MOI** — the Poisson distribution models random lentiviral integration events, and the zero-truncated correction accounts for only observing cells that received at least one barcode
- **The estimated λ should approximate your MOI** — if you infected at MOI ~6, you should see λ ≈ 6 in the fit results
- **A poor fit may indicate**: doublets, barcode cross-talk, insufficient noise filtering, or non-random integration

### Potential limitations and considerations

1. **Noise sensitivity**: At high MOI, each cell has more barcode entries, increasing the chance that noise barcodes survive filtering. Consider using a **stricter Z-threshold** (2.5–3.0) and inspecting the filtered output carefully.

2. **Clone resolution**: With 6–7 barcodes per cell, clone barcodes become long concatenated strings. Two cells from the same clone must share the **exact same set** of barcodes to be grouped. If noise filtering removes a real barcode from some cells but not others, clonal cells may receive different clone IDs. This is inherently harder to resolve at high MOI.

3. **Starcode clustering (Step 8)**: At high MOI, cells have more lineage barcode sequences to cluster. If different real barcodes are within the Starcode distance threshold of each other, they may be incorrectly merged. Consider **reducing `--max_distance`** (e.g., from 3 to 2) if your barcode library has limited diversity.

4. **`--max-lineage-count`**: Strongly recommended for high MOI experiments. Cells with extremely high lineage counts are likely doublets or other artifacts. Set this to 2–3× your expected MOI to exclude them from clone assignment while retaining real multi-integration cells.

5. **Computational cost**: Starcode (Step 8) runs per-cell, so cells with many barcodes generate larger inputs. For very high MOI with large cell counts, this step may take longer.

---

## Notes

- The pipeline is designed for **10x Genomics Chromium** (v2/v3) and **Drop-seq** based REWIND experiments.
- Steps 4 and 5 (`collapse.py` and `addressct2.py`) handle a specific artifact of oligo-dT primed barcoding where incomplete extension causes barcode length errors — this may not be needed for all experimental setups.
- The final clone assignment step assumes most cells carry a **single lentiviral integration** (`--top-n 1`) by default, fitting a truncated Poisson model to validate this assumption. See the [High MOI Scenarios](#high-moi-multiple-integrations-per-cell-scenarios) section for multi-integration experiments.
- Output file names are derived from input file paths in several scripts; organize your working directory accordingly.
- The `run_pipeline.py` master script handles file naming automatically when used.

---

## License

Please refer to the repository owner for licensing information.
