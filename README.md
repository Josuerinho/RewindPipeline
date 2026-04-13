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
**Purpose:** Extracts cell barcodes (16 nt), UMIs (12 nt), and lineage barcodes from paired-end FASTQ files (10x Genomics format). Locates lineage barcodes in Read 2 by searching for a configurable capture sequence (default: `GCTCACCTATTAGCGGCTAAGG`) and extracting a configurable number of nucleotides upstream (default: 84 nt). Counts mismatches to the expected **WSN trinucleotide repeat pattern** (W = A/T, S = G/C, N = any).

```bash
python rewind_extractor.py \
  --fastq1 <read1.fastq.gz> \
  --fastq2 <read2.fastq.gz> \
  --output <rewind_address.tsv> \
  --instrument-run-flowcell-ID <e.g., AV100007:PY055:2336402118:> \
  --capture-seq GCTCACCTATTAGCGGCTAAGG \
  --barcode-length 84
```

#### Understanding `--instrument-run-flowcell-ID`

This parameter is the **prefix string** from your FASTQ read headers that identifies the instrument, run, and flowcell. It is used as a **string delimiter** in `record.id.split(prefix)` to strip the prefix and extract the unique lane:tile:x:y portion of each read ID.

A typical FASTQ read ID looks like:

```
@AV100007:PY055:2336402118:1:1101:3425:1000
```

This has the structure: `instrument:run:flowcell:lane:tile:x:y`

By passing `--instrument-run-flowcell-ID "AV100007:PY055:2336402118:"` (note the trailing colon), the code performs:

```python
record.id.split("AV100007:PY055:2336402118:")
# → ["", "1:1101:3425:1000"]
```

This extracts `"1:1101:3425:1000"` (the lane:tile:x:y portion) as the unique read identifier. The trailing colon ensures a clean split. You should copy this prefix directly from the first line of your FASTQ file, including the trailing colon.

#### Configurable capture sequence and barcode length

- **`--capture-seq`**: The DNA sequence in Read 2 used to locate the lineage barcode. The barcode is extracted immediately upstream of this sequence. Change this if using a different plasmid backbone. Default: `GCTCACCTATTAGCGGCTAAGG`.
- **`--barcode-length`**: Number of nucleotides upstream of the capture sequence to extract as the lineage barcode. Increase this if your barcode construct is longer than 84 nt. Default: `84`.

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
| **Instrument/Run/Flowcell ID** | Header prefix from FASTQ read IDs used to extract lane:tile:x:y (e.g., `AV100007:PY055:2336402118:`) — see [Step 1 docs](#understanding---instrument-run-flowcell-id) for details |
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
| `--capture-seq` | `rewind_extractor.py` | `GCTCACCTATTAGCGGCTAAGG` | Capture sequence for locating lineage barcodes in Read 2 |
| `--barcode-length` | `rewind_extractor.py` | `84` | Lineage barcode length (nt) upstream of capture sequence |
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

### Custom capture sequence and barcode length

If your plasmid uses a different capture sequence or a longer barcode region:

```bash
python run_pipeline.py \
  --fastq1 sample_R1.fastq.gz \
  --fastq2 sample_R2.fastq.gz \
  --instrument-run-flowcell-ID "AV100007:PY055:2336402118:" \
  --technology 10xv3 \
  --whitelist-file barcode_translation.tsv.gz \
  --edrops-file edrops_filtered_cells.tsv \
  --starcode-path /path/to/starcode \
  --capture-seq ATCGATCGATCGATCGATCG \
  --barcode-length 120
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
| `--capture-seq` | `GCTCACCTATTAGCGGCTAAGG` | Capture sequence for locating lineage barcodes in Read 2 |
| `--barcode-length` | `84` | Lineage barcode length (nt) upstream of capture sequence |
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

## Timepoint Matching Analysis

When you have data from two experimental timepoints (e.g., T0 baseline and T1 after treatment), you can match T1 clones back to their T0 ancestors using the **`match_timepoints.py`** script. This is useful for tracking clonal dynamics such as drug resistance or differential expansion.

### Algorithm

1. **Load** clone → barcode mappings from both timepoints (`*_cloneID_cloneBarcode.tsv` from Step 10).
2. **Fuzzy barcode matching**: For each barcode in a T1 clone, find all T0 barcodes within a configurable Hamming distance threshold (default: 2). This accounts for sequencing errors and low-level barcode mutations.
3. **Candidate scoring**: For each (T1 clone, T0 clone) candidate pair, compute a **Fuzzy Jaccard similarity**:
   ```
   similarity = |matched_barcodes| / |T1_barcodes ∪ T0_barcodes|
   ```
   Optionally, similarity can be weighted by UMI counts to prioritise high-confidence barcodes.
4. **Ranking and filtering**: Rank T0 candidates per T1 clone; report top-N matches and flag T1 clones below the similarity threshold as unmatched.

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--hamming-threshold` | 2 | Max Hamming distance for fuzzy barcode matching |
| `--similarity-threshold` | 0.8 | Min Fuzzy Jaccard score to accept a match |
| `--top-n-matches` | 5 | Number of ranked T0 candidates per T1 clone |
| `--weight-by-umi` | off | Weight similarity by UMI counts |
| `--n-jobs` | 1 | Parallel workers (0 = auto-detect) |

### Standalone Usage

```bash
# Basic matching
python match_timepoints.py \
    --t0-clone-barcode T0_cloneID_cloneBarcode.tsv \
    --t1-clone-barcode T1_cloneID_cloneBarcode.tsv \
    --output-dir timepoint_matching_results

# With UMI weighting and relaxed Hamming threshold
python match_timepoints.py \
    --t0-clone-barcode T0_cloneID_cloneBarcode.tsv \
    --t1-clone-barcode T1_cloneID_cloneBarcode.tsv \
    --t0-filtered T0_rewind_filtered.tsv \
    --t1-filtered T1_rewind_filtered.tsv \
    --weight-by-umi \
    --hamming-threshold 3 \
    --similarity-threshold 0.85 \
    --n-jobs 4 \
    --output-dir timepoint_matching_results
```

### Integrated Pipeline Usage

Add `--run-timepoint-matching` and `--t0-reference` to `run_pipeline.py`:

```bash
python run_pipeline.py \
    --fastq1 T1_R1.fastq.gz \
    --fastq2 T1_R2.fastq.gz \
    --instrument-run-flowcell-ID "INSTRUMENT:RUN:FLOWCELL:" \
    --technology 10xv3 \
    --whitelist-file barcode_translation.tsv.gz \
    --edrops-file edrops_filtered_cells.tsv \
    --starcode-path /path/to/starcode \
    --output-dir T1_results/ \
    --sample-name T1_sample \
    --run-timepoint-matching \
    --t0-reference /path/to/T0_cloneID_cloneBarcode.tsv \
    --stop-step 11
```

### Output Files

| File | Description |
|---|---|
| `best_matches.txt` | Best T0 match per T1 clone (T1_clone_id, T0_clone_id, similarity, counts) |
| `ranked_matches.txt` | Top-N T0 matches per T1 clone with scores |
| `unmatched_clones.txt` | T1 clones whose best match is below the similarity threshold |
| `matching_statistics.txt` | Summary: match rate, mean/median/std similarity, etc. |
| `plots/` | QC visualisations (see below) |

### QC Plots

- **similarity_histogram.png** — Distribution of best-match scores with threshold line
- **clone_size_scatter.png** — T0 vs T1 clone size coloured by similarity
- **matched_vs_unmatched.png** — Bar chart of match status counts
- **top_matches_heatmap.png** — Heatmap of highest-scoring clone pairs
- **similarity_vs_matched_count.png** — Scatter of similarity vs. number of shared barcodes

### Interpreting Results

- **High similarity (≥ 0.9)**: Confident one-to-one match. The T1 clone is almost certainly descended from that T0 clone.
- **Moderate similarity (0.7–0.9)**: Likely match but with barcode dropout. Inspect the number of matched barcodes — higher counts increase confidence.
- **Low similarity (< 0.7)**: Ambiguous. Multiple T0 candidates may share a few barcodes. Consider these matches provisional.
- **Unmatched clones**: May indicate very high dropout, novel barcode combinations from recombination, or rare clones whose T0 ancestor had too few cells to be detected.
- **Multiple high-scoring T0 matches**: Could indicate closely related clones in T0 or barcode sharing artifacts. The ranked output helps assess this.

### Scalability

The algorithm is designed for large datasets (80K × 20K comparisons):
- Uses dictionary-based indexing to avoid exhaustive pairwise comparisons
- Only computes similarity for T0 clones that share at least one fuzzy-matched barcode with the T1 clone
- Supports multiprocessing via `--n-jobs` for CPU-parallel matching
- Progress bars report completion during long runs

---

## Differential Enrichment Analysis

### Overview

The **Differential Enrichment Analysis** (Step 12) identifies clones that are significantly enriched or depleted due to drug selection by comparing the Treated/T0 fold-change against the Control/T0 fold-change. Using the untreated Control population as a reference accounts for **random clonal drift** (stochastic changes in clone frequency during the 2-month passage), ensuring that only drug-specific effects are reported as significant.

### Statistical Framework

For every clone observed at T0, the analysis computes:

1. **Fold-changes** (with a configurable pseudocount, default = 1):
   - `FC_control = (Control_size + pseudocount) / (T0_size + pseudocount)`
   - `FC_treated = (Treated_size + pseudocount) / (T0_size + pseudocount)`

2. **Log2 Fold-Change Ratio** — the key metric:
   - `Log2_FC_ratio = log2(FC_treated / FC_control)`
   - A positive value means the clone grew more under drug treatment than in the control.
   - A negative value means the clone shrank more under drug treatment.

3. **Fisher's exact test** for each clone on the 2×2 contingency table `[[T0+1, Control+1], [T0+1, Treated+1]]` to assess whether the difference in fold-changes is statistically significant.

4. **Benjamini-Hochberg FDR correction** to control false discoveries across many simultaneous tests.

5. **Classification** (configurable thresholds):
   - **Enriched**: q < 0.05 AND Log2_FC_ratio > 1.0
   - **Depleted**: q < 0.05 AND Log2_FC_ratio < −1.0
   - **Neutral**: everything else

### Why the Control Is Essential

Without a Control population, any observed change in clone frequency from T0→T1 could be due to:
- Random drift during passage
- Differences in sequencing depth
- Sampling noise

By normalising Treated fold-changes against Control fold-changes, these confounders cancel out and the remaining signal reflects **drug-specific** selection.

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pseudocount` | 1 | Added to numerator and denominator before fold-change calculation. Prevents division by zero and stabilises ratios for small clones. |
| `--fdr-threshold` | 0.05 | q-value cutoff for significance. |
| `--log2fc-threshold` | 1.0 | Minimum absolute Log2 FC ratio required (in addition to FDR) to call a clone enriched/depleted. A value of 1.0 means the clone must show at least a 2-fold difference between Treated and Control. |

### Standalone Usage

After running the pipeline (Steps 1–10) on all three populations (T0, Control, Treated) and performing timepoint matching (Step 11), run enrichment analysis independently:

```bash
python differential_enrichment_analysis.py \
    --t0-clones T0_output/rewind_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --control-clones Control_output/rewind_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --treated-clones Treated_output/rewind_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --control-matches control_matching/best_matches.txt \
    --treated-matches treated_matching/best_matches.txt \
    --output-dir enrichment_results
```

With UMI-based clone sizing (recommended when rewind_filtered.tsv files are available):

```bash
python differential_enrichment_analysis.py \
    --t0-clones T0_cloneID_cloneBarcode.tsv \
    --control-clones Control_cloneID_cloneBarcode.tsv \
    --treated-clones Treated_cloneID_cloneBarcode.tsv \
    --control-matches control_matching/best_matches.txt \
    --treated-matches treated_matching/best_matches.txt \
    --t0-filtered T0_rewind_filtered.tsv \
    --control-filtered Control_rewind_filtered.tsv \
    --treated-filtered Treated_rewind_filtered.tsv \
    --output-dir enrichment_results \
    --pseudocount 1 --fdr-threshold 0.05 --log2fc-threshold 1.0
```

### Integrated Pipeline Usage

Run enrichment analysis as part of the pipeline (Steps 1–12). This is useful when the current pipeline run is the Treated sample:

```bash
python run_pipeline.py \
    --fastq1 treated_R1.fastq.gz --fastq2 treated_R2.fastq.gz \
    --instrument-run-flowcell-ID "INSTRUMENT:RUN:FLOWCELL:" \
    --technology 10xv3 \
    --whitelist-file barcode_translation.tsv.gz \
    --edrops-file edrops_filtered_cells.tsv \
    --starcode-path /path/to/starcode \
    --output-dir treated_output \
    --run-timepoint-matching --t0-reference T0_cloneID_cloneBarcode.tsv \
    --run-enrichment-analysis \
    --control-reference Control_cloneID_cloneBarcode.tsv \
    --control-matches control_matching/best_matches.txt \
    --stop-step 12
```

### Output Files

| File | Description |
|------|-------------|
| `differential_enrichment_results.txt` | Full results table with Clone_ID, T0_size, Control_size, Treated_size, FC_control, FC_treated, Log2_FC_ratio, p_value, q_value, Status |
| `enriched_clones.txt` | Subset: significantly enriched clones only |
| `depleted_clones.txt` | Subset: significantly depleted clones only |
| `neutral_clones.txt` | Subset: non-significant clones |
| `summary_statistics.txt` | Overview statistics and top enriched/depleted clones |
| `plots/` | Directory with publication-quality visualisation plots |

### QC Plots

1. **Volcano plot** (`volcano_plot.png`): Log2_FC_ratio (x) vs −log10(q-value) (y). Enriched clones appear in the upper-right (red), depleted in the upper-left (blue). Dashed lines mark significance thresholds.

2. **MA plot** (`ma_plot.png`): Log2 mean abundance (x) vs Log2_FC_ratio (y). Shows whether enrichment is related to clone size — important for spotting artefacts in very small or very large clones.

3. **Clone trajectory plot** (`clone_trajectory_plot.png`): Log2 FC Control (x) vs Log2 FC Treated (y). Points on the diagonal grew equally in both conditions; points above the diagonal are preferentially enriched under treatment.

4. **Enrichment distribution** (`enrichment_distribution.png`): Histogram of all Log2_FC_ratio values. A symmetric distribution centred near zero suggests most clones behave similarly; heavy tails indicate strong selection.

5. **Top clones heatmap** (`top_clones_heatmap.png`): Log2 UMI counts for the top 20 enriched and top 20 depleted clones across T0, Control, and Treated.

6. **Summary bar plot** (`summary_bar_plot.png`): Counts of Enriched / Neutral / Depleted clones.

### Interpreting Results

#### How to read the Volcano plot
- **Upper-right quadrant** (high −log10(q), high Log2_FC): Clones that are both statistically significant and biologically meaningful — **true resistance candidates**.
- **Upper-left quadrant**: Clones significantly **depleted** by the drug.
- **Bottom band** (low −log10(q)): Clones without enough evidence for differential behaviour.

#### What Log2_FC_ratio values mean
| Log2_FC_ratio | Interpretation |
|---------------|----------------|
| > 2 | Clone grew >4× more in Treated vs Control — strong enrichment |
| 1 to 2 | 2–4× enrichment — moderate |
| −1 to 1 | Less than 2-fold difference — likely noise or drift |
| −2 to −1 | 2–4× depletion — moderate sensitivity |
| < −2 | >4× depletion — strong drug sensitivity |

#### Identifying true resistance clones
1. Start with the **enriched_clones.txt** file.
2. Prioritise clones with high Log2_FC_ratio **and** low q-value.
3. Cross-reference with the MA plot — very small clones with extreme fold-changes may be unreliable.
4. Check the clone trajectory plot to confirm the clone grew specifically in the Treated arm.

#### Common pitfalls
- **Too many significant results?** Tighten thresholds: increase `--log2fc-threshold` or decrease `--fdr-threshold`.
- **No significant results?** Check that timepoint matching worked well (inspect `matching_statistics.txt`). Consider relaxing thresholds or verifying that the drug is effective.
- **Batch effects**: Ensure T0, Control, and Treated were processed with the same pipeline parameters.
- **Very small clones**: A single-cell clone that drops out can produce extreme fold-changes; interpret with caution. The pseudocount helps, but small clones are inherently noisy.

### Parameter Tuning Recommendations

- **pseudocount**: Increase to 5–10 if you have many small clones and want to reduce extreme fold-changes in low-abundance clones. Keep at 1 for standard analysis.
- **fdr-threshold**: 0.05 is standard. Use 0.01 for a more conservative analysis.
- **log2fc-threshold**: 1.0 (2-fold) is a common biological cutoff. Use 0.585 (~1.5-fold) for more permissive analysis or 2.0 (4-fold) for strict filtering.

### Example Workflow: FASTQ to Enrichment

```bash
# 1. Process T0 sample (Steps 1–10)
python run_pipeline.py --fastq1 T0_R1.fq.gz --fastq2 T0_R2.fq.gz \
    [pipeline params] --output-dir T0_output --sample-name T0

# 2. Process Control sample (Steps 1–10)
python run_pipeline.py --fastq1 Ctrl_R1.fq.gz --fastq2 Ctrl_R2.fq.gz \
    [pipeline params] --output-dir Ctrl_output --sample-name Ctrl

# 3. Process Treated sample (Steps 1–10)
python run_pipeline.py --fastq1 Treat_R1.fq.gz --fastq2 Treat_R2.fq.gz \
    [pipeline params] --output-dir Treat_output --sample-name Treat

# 4. Run timepoint matching for Control → T0
python match_timepoints.py \
    --t0-clone-barcode T0_output/T0_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --t1-clone-barcode Ctrl_output/Ctrl_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --output-dir ctrl_matching

# 5. Run timepoint matching for Treated → T0
python match_timepoints.py \
    --t0-clone-barcode T0_output/T0_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --t1-clone-barcode Treat_output/Treat_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --output-dir treat_matching

# 6. Run differential enrichment analysis
python differential_enrichment_analysis.py \
    --t0-clones T0_output/T0_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --control-clones Ctrl_output/Ctrl_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --treated-clones Treat_output/Treat_corrected_UMI_counts.z2.cloneID_cloneBarcode.tsv \
    --control-matches ctrl_matching/best_matches.txt \
    --treated-matches treat_matching/best_matches.txt \
    --t0-filtered T0_output/T0_corrected_UMI_counts.z2.rewind_filtered.tsv \
    --control-filtered Ctrl_output/Ctrl_corrected_UMI_counts.z2.rewind_filtered.tsv \
    --treated-filtered Treat_output/Treat_corrected_UMI_counts.z2.rewind_filtered.tsv \
    --output-dir enrichment_results
```

---

## License

Please refer to the repository owner for licensing information.
