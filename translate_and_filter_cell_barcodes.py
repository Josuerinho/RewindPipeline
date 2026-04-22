#!/usr/bin/env python3

import pandas as pd
import argparse

def cellbcs_count(df, cellbc_column = "cell_barcode_gex"):
    return len(df[cellbc_column].unique())

def read_rewind_file(file_path):
    """
    Read the rewind_address.max_mismatches0.fcts.txt file into a DataFrame.

    Parameters
    ----------
    file_path : str
        Path to the rewind_address.max_mismatches0.fcts.txt file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the data from the file.
    """
    return pd.read_csv(file_path, sep='\t', header=None, names=['cell_barcode_feature', 'umi', 'lineage_barcode', 'read_counts'])

def read_whitelist(file_path):
    """
    Read the whitelist file into a DataFrame.

    Parameters
    ----------
    file_path : str
        Path to the whitelist file; delimiter = "\t"

    Returns
    -------
    pd.DataFrame
        DataFrame containing the whitelist data.
    """
    file_extension = file_path.rsplit(".", 1)[-1]
    if file_extension == "gz":
        return pd.read_csv(file_path, compression = "gzip", delimiter = "\t",  header=None, names=['cell_barcode_gex', 'cell_barcode_feature'])
    else:
        return pd.read_csv(file_path, sep='\t', header=None, names=['cell_barcode_gex', 'cell_barcode_feature'])

def translate_barcodes(rewind_df, whitelist_df):
    """
    Translate cell barcodes in the rewind DataFrame using the whitelist DataFrame.

    Parameters
    ----------
    rewind_df : pd.DataFrame
        DataFrame containing the rewind data.
    whitelist_df : pd.DataFrame
        DataFrame containing the whitelist data.

    Returns
    -------
    pd.DataFrame
        DataFrame with translated cell barcodes.
    """
    print(f'number of cell_barcode_feature: {cellbcs_count(rewind_df, "cell_barcode_feature")}')

    translation_dict = dict(zip(whitelist_df['cell_barcode_feature'], whitelist_df['cell_barcode_gex']))
    cell_barcode_gex = rewind_df['cell_barcode_feature'].map(translation_dict)
    rewind_df.insert(0, "cell_barcode_gex", cell_barcode_gex)
    rewind_df.drop(columns = "cell_barcode_feature", inplace = True)
    rewind_df = rewind_df.dropna().reset_index(drop = True)

    print(f'number of cell_barcode_gex: {cellbcs_count(rewind_df)}')
    return rewind_df

def filter_barcodes(rewind_df, edrops_file_path):
    """
    Filter translated cell barcodes using edrops filtered cell barcodes.

    Parameters
    ----------
    rewind_df : pd.DataFrame
        DataFrame with translated cell barcodes.
    edrops_file_path : str
        Path to the edrops filtered cell barcodes file. It can also be just the cell barcodes in a single column with no header.

    Returns
    -------
    pd.DataFrame
        DataFrame with filtered cell barcodes.
    """
    
    # adjust to accept the edrops file with just the cell barcodes in a single column with no header
    edrops_df = pd.read_csv(edrops_file_path, sep='\t', header=None)
    if edrops_df.shape[1] == 1:
        edrops_df.columns = ['cell_barcode_gex']
    elif edrops_df.shape[1] == 3:
        edrops_df.columns = ['cell_barcode_gex', "molecules_count", "genes_count"]
    else:
        raise ValueError("edrops file should have either 1 column (cell_barcode_gex) or 3 columns (cell_barcode_gex, molecules_count, genes_count)")
    print(f'number of cell_barcode_gex in edrops: {edrops_df.shape[0]}')

    filtered_rewind_df = rewind_df[rewind_df['cell_barcode_gex'].isin(edrops_df['cell_barcode_gex'])]
    filtered_rewind_df = filtered_rewind_df.reset_index(drop = True)
    print(f'number of edrops filtered cell_barcode_gex: {cellbcs_count(filtered_rewind_df)}')
    print(f'{cellbcs_count(filtered_rewind_df) / edrops_df.shape[0]}')

    return filtered_rewind_df

def main(rewind_file_path, whitelist_file_path, edrops_file_path, output_file_path, min_reads_count):
    """
    Main function to translate and filter cell barcodes.

    Parameters
    ----------
    rewind_file_path : str
        Path to the rewind_address.max_mismatches0.fcts.txt file.
    whitelist_file_path : str
        Path to the whitelist file.
    edrops_file_path : str
        Path to the edrops filtered cell barcodes file.
    output_file_path : str
        Path to the output file where the filtered data will be saved.
    """
    rewind_df = read_rewind_file(rewind_file_path)
    rewind_df = rewind_df[rewind_df['read_counts'] > min_reads_count]
    whitelist_df = read_whitelist(whitelist_file_path)
    translated_rewind_df = translate_barcodes(rewind_df, whitelist_df)
    filtered_rewind_df = filter_barcodes(translated_rewind_df, edrops_file_path)
    filtered_rewind_df.to_csv(output_file_path, sep='\t', index=False, header = None)
    print(f"Filtered data saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate and filter cell barcodes.")
    parser.add_argument("--rewind_file", required=True, help="Path to the rewind_address.max_mismatches0.fcts.txt file.")
    parser.add_argument("--whitelist_file", required=True, help="Path to the whitelist file.")
    parser.add_argument("--edrops_file", required=True, help="Path to the edrops filtered cell barcodes file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file where the filtered data will be saved.")
    parser.add_argument("--min-reads-count", required=True, type=int, help="Molecules need to have greater than this minimum of reads to be included for downstream process") 
    args = parser.parse_args()
    main(args.rewind_file, args.whitelist_file, args.edrops_file, args.output_file, args.min_reads_count)

