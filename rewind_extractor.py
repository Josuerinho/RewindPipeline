#!/usr/bin/env python3
import argparse
import gzip
from Bio import SeqIO

DEFAULT_CAPTURE_SEQ = "GCTCACCTATTAGCGGCTAAGG"
DEFAULT_BARCODE_LENGTH = 84

# Define the sets for W, S, and N conditions
W = set('AT')
S = set('GC')
N = set('ATGC')  # Though we won't actually check N as all are allowed

# Function to compute mismatches against the expected WSN pattern
def count_wsn_mismatches(sequence):
    mismatches = 0
    for i in range(0, len(sequence), 3):
        if i < len(sequence) and sequence[i] not in W:
            mismatches += 1
        if i + 1 < len(sequence) and sequence[i + 1] not in S:
            mismatches += 1
        # N is skipped as it can be any nucleotide
    return mismatches

def rewind_extractor(fastq1_path, fastq2_path, output_path, instrument_run_flowcell_ID,
                     capture_seq=DEFAULT_CAPTURE_SEQ, barcode_length=DEFAULT_BARCODE_LENGTH):
    """
    Extracts cell barcodes, UMIs, and lineage barcodes from read1 and read2 FASTQ files, respectively,
    and writes them to an output file along with the number of mismatches to the WSN pattern.

    Parameters
    ----------
    fastq1_path : str
        Path to the read1 FASTQ file which contains cell barcodes and UMIs.
    fastq2_path : str
        Path to the read2 FASTQ file which contains the lineage barcodes.
    output_path : str
        Path to the output file where results will be saved.
    instrument_run_flowcell_ID : str
        The instrument:run:flowcell prefix string from the FASTQ read header, used to
        split record.id and extract the lane:tile:x:y read identifier portion.
    capture_seq : str
        The capture sequence used to locate lineage barcodes in Read 2.
        Default: 'GCTCACCTATTAGCGGCTAAGG'
    barcode_length : int
        Length of the lineage barcode (in nt) upstream of the capture sequence.
        Default: 84

    Output
    ------
    The output file will contain the following tab-separated fields:
    - readid: The identifier for the read.
    - cell_barcode: The cell barcode extracted from read1.
    - umi: The UMI extracted from read1.
    - lineage_barcode: The lineage barcode extracted from read2.
    - mismatches: The number of mismatches between the lineage barcode and the expected WSN pattern.

    Example
    -------
    Given an input read1 FASTQ file and read2 FASTQ file, the function extracts relevant data and writes it
    to the output file with the following format:
    readid    cell_barcode    umi    lineage_barcode    mismatches
    """
    cell_barcode_len = 16
    umi_len = 12

    # Lists to store cell barcodes, UMIs, and read IDs
    cell_barcodes = []
    umis = []
    read_ids = []

    print("Extract from Read1 ------------------------------")
    # Extract cell barcodes, UMIs, and read IDs from read1
    with gzip.open(fastq1_path, 'rt') as f1:
        for record in SeqIO.parse(f1, "fastq"):
            read_id = record.id.split(instrument_run_flowcell_ID)[1]
            cell_barcode = str(record.seq[:cell_barcode_len])
            umi = str(record.seq[cell_barcode_len:cell_barcode_len + umi_len])
            cell_barcodes.append(cell_barcode)
            umis.append(umi)
            read_ids.append(read_id)

    print("Extract from Read2 ------------------------------")
    # Extract lineage barcodes from read2 and write to the output file
    number_of_readids_unmatched = 0
    number_of_lineage_barcodes_mismatched = 0
    number_of_reads_with_captureseq2 = 0
    number_reads_in_address = 0
    with gzip.open(fastq2_path, 'rt') as f2, open(output_path, 'w') as out:
        for i, record in enumerate(SeqIO.parse(f2, "fastq")):
            read_id = record.id.split(instrument_run_flowcell_ID)[1]
            read2_sequence = str(record.seq)
            capture_seq2_index = read2_sequence.find(capture_seq)
            if capture_seq2_index != -1:
                number_of_reads_with_captureseq2 += 1
                lineage_barcode_start = capture_seq2_index - barcode_length
                if lineage_barcode_start >= 0:
                    lineage_barcode = read2_sequence[lineage_barcode_start:capture_seq2_index]
                    mismatches = count_wsn_mismatches(lineage_barcode)
                    if mismatches > 0:
                        number_of_lineage_barcodes_mismatched += 1
                    if read_id == read_ids[i]:
                        out.write(f"{read_id}\t{cell_barcodes[i]}\t{umis[i]}\t{lineage_barcode}\t{mismatches}\n")
                        number_reads_in_address += 1
                    else:
                        number_of_readids_unmatched += 1

    print(f'Number of lineage barcodes with at least one mismatch to the WSN pattern: {number_of_lineage_barcodes_mismatched}')
    print(f'Number of reads in Address: {number_reads_in_address}')
    print(f'Fraction of mismatched lineage barcode: {number_of_lineage_barcodes_mismatched/number_reads_in_address}')
    print(f'Number of reads containing CaptureSeq2: {number_of_reads_with_captureseq2}')
    print(f'Number of readids did not match between Read1 and Read2: {number_of_readids_unmatched}')
    print(f"Extraction complete. Results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Extract cell barcodes, UMIs, and lineage barcodes from FASTQ files.')
    parser.add_argument('--fastq1', type=str, required=True, help='Path to the read1 FASTQ file')
    parser.add_argument('--fastq2', type=str, required=True, help='Path to the read2 FASTQ file')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file')
    parser.add_argument('--instrument-run-flowcell-ID', type=str, required=True,
                        help='Instrument:Run:Flowcell ID prefix from the FASTQ read header, '
                             'used to split read IDs. Example: "AV100007:PY055:2336402118:" — '
                             'the code calls record.id.split(this_value) to extract the '
                             'lane:tile:x:y portion after this prefix.')
    parser.add_argument('--capture-seq', type=str, default=DEFAULT_CAPTURE_SEQ,
                        help='Capture sequence used to locate lineage barcodes in Read 2. '
                             f'(default: {DEFAULT_CAPTURE_SEQ})')
    parser.add_argument('--barcode-length', type=int, default=DEFAULT_BARCODE_LENGTH,
                        help='Length of the lineage barcode (nt) upstream of the capture sequence. '
                             f'(default: {DEFAULT_BARCODE_LENGTH})')
    args = parser.parse_args()

    rewind_extractor(args.fastq1, args.fastq2, args.output, args.instrument_run_flowcell_ID,
                     capture_seq=args.capture_seq, barcode_length=args.barcode_length)

if __name__ == "__main__":
    main()

