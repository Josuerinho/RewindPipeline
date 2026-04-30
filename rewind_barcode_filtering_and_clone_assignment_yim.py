import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.optimize import fsolve
from scipy.stats import poisson, kstest
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os


def parse_args():
    """
    Parses command-line arguments for the rewind_barcode_filtering_and_cloneid_assignment Script.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="This script processes lineage barcode data to filter out noise, analyze clone barcodes, and produce various QC plots. The input is a TSV file containing cell and lineage barcodes. The script outputs filtered data, statistical summaries, and QC plots.")
    parser.add_argument( '--input-file', type=str, required=True, help="Path to the input TSV file containing starcode corrected rewind lineage barcode address data. This file should have columns for cell_barcode and corrected_lineage_barcode. The input path before the file extension would be used as the output directory path and prefix.")
    parser.add_argument( '--top-n', type=int, default=1, help="Number of top log-normalized counts to consider for calculating the expected count for each cell. Default is 1, because assuming one lineage barcode per cell so the top molecules count lineag barcode should be real.")
    parser.add_argument( '--z-threshold', type=float, default=2, help="Z-score threshold for filtering out noise. Rows with Z-scores above this threshold will be considered noise. Default is 2.")
    parser.add_argument( '--min-group-size', type=int, default=5, help="Minimum group size to avoid handling small groups separately. Groups smaller than this size will be flagged. Default is 5.")
    parser.add_argument( '--max-lineage-count', type=int, default=None, help="Maximum lineage count to consider for determining clone barcodes. Cells with lineage counts above this value will be excluded from clone barcode analysis. If not specified, all lineage counts will be considered.")
    return parser.parse_args()


def filter_noise(df, top_n, z_threshold, min_group_size, target_sum=1e4):
    """
    Filters out noise in lineage barcodes molecule counts by normalizing, log-transforming,
    and calculating Z-scores to identify significant deviations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the lineage barcode data with columns 'cell_barcode', 'umi_count', and 'corrected_lineage_barcode'.
    top_n : int
        Number of top log-normalized counts to consider for calculating the expected count for each cell.
    z_threshold : float
        Z-score threshold for filtering out noise. Rows with Z-scores above this threshold will be considered noise.
    min_group_size : int
        Minimum group size to avoid handling small groups separately. Groups smaller than this size will be flagged.
    target_sum : float, optional
        Target sum for normalization of UMI counts. Default is 1e4.

    Returns
    -------
    tuple
        A tuple containing:
        - df (pd.DataFrame): Original DataFrame with added normalization, log transformation, and Z-score columns.
        - filtered_df (pd.DataFrame): DataFrame with noise filtered out based on the Z-score threshold.
    """

    # Normalize UMI counts for each cell to sum to target_sum
    df['normalized_umi_count'] = df.groupby('cell_barcode')['umi_count'].transform(lambda x: (x / x.sum()) * target_sum)

    # Log-transform the normalized counts
    df['log_normalized_umi_count'] = np.log1p(df['normalized_umi_count'])

    # Calculate the mean of the top N highest log-normalized counts as the expected count for each cell
    df['rank'] = df.groupby('cell_barcode')['log_normalized_umi_count'].rank(method='first', ascending=False)
    expected_counts = df[df['rank'] <= top_n].groupby('cell_barcode')['log_normalized_umi_count'].mean().reset_index()
    expected_counts = expected_counts.rename(columns={'log_normalized_umi_count': 'expected_count'})

    # Merge expected counts back to the main DataFrame
    df = pd.merge(df, expected_counts, on='cell_barcode')

    # Compute the standard deviation for cells with the same number of unique lineage barcodes: subpopulation standard deviation
    lineage_count_per_cell = df.groupby('cell_barcode')['corrected_lineage_barcode'].nunique().reset_index(name='lineage_count')
    df = pd.merge(df, lineage_count_per_cell, on='cell_barcode')
    std_devs = df.groupby('lineage_count')['log_normalized_umi_count'].std().reset_index(name='std_dev')

    # Merge standard deviations back to the main DataFrame
    df = pd.merge(df, std_devs, on='lineage_count', how='left')

    # Identify small subpopulations
    group_sizes = df.groupby('cell_barcode')['lineage_count'].transform('size')
    df['is_small_group'] = group_sizes < min_group_size

    # Compute the Z-score for log-normalized counts relative to the expected count and subpopulation std dev
    def calculate_z_score(row):
        if row['is_small_group']:
            return np.nan  # Indicating insufficient data due to small group size
        elif row['std_dev'] == 0 and row['lineage_count'] == 1:
            return 0  # When std_dev is zero and lineage_count is 1
        else:
            return (row['log_normalized_umi_count'] - row['expected_count']) / row['std_dev']

    df['z_score'] = df.apply(calculate_z_score, axis=1)

    # Filter out noise based on Z-score
    filtered_df = df[~df['z_score'].isna() & (df['z_score'].abs() < z_threshold)]

    # Update lineage_count in the filtered DataFrame
    filtered_df['lineage_count'] = filtered_df.groupby('cell_barcode')['corrected_lineage_barcode'].transform('nunique')

    return df, filtered_df

def test_truncated_poisson_distribution(total_df, column_name):
    """
    Test if the distribution of total lineage barcodes per cell fits a truncated Poisson distribution, where P(X=0) is removed

    Parameters
    ----------
    total_df : pd.DataFrame
        DataFrame containing the cell_barcode and lineages_per_cell.
    column_name : str
        Name of the column containing the total number of lineage barcodes per cell.

    Returns
    -------
    tuple
        A tuple containing:
        - empirical_counts (np.ndarray): Empirical counts of lineage barcodes per cell.
        - expected_counts (np.ndarray): Expected counts from the truncated Poisson distribution.
        - expected_cdf (np.ndarray): Cumulative distribution function (CDF) of the expected counts.
        - empirical_cdf (np.ndarray): CDF of the empirical counts.
        - results (dict): Dictionary containing observed mean, estimated lambda, K-S statistic, p-value, and percentage of cells with at least one lentivirus.
    """

    # Step 1: Calculate the empirical distribution
    empirical_counts = total_df[column_name].value_counts().sort_index()
    empirical_cdf = np.cumsum(empirical_counts / empirical_counts.sum())

    # Step 2: Calculate the observed mean
    observed_mean = total_df[column_name].mean()

    # Step 3: Define the function to solve for λ (lambda)
    def solve_lambda(lambda_estimate):
        return lambda_estimate / (1 - np.exp(-lambda_estimate)) - observed_mean

    # Use fsolve to find the root of the equation for λ
    lambda_estimate_initial = observed_mean  # Initial guess
    lambda_estimate = fsolve(solve_lambda, lambda_estimate_initial)[0]

    # Step 4: Adjust for the truncated Poisson distribution
    truncation_correction = 1 - poisson.pmf(0, lambda_estimate)
    max_lineages = empirical_counts.index.max()
    poisson_dist = poisson.pmf(k=np.arange(1, max_lineages + 1), mu=lambda_estimate) / truncation_correction
    expected_counts = poisson_dist * empirical_counts.sum()
    expected_cdf = np.cumsum(expected_counts) / empirical_counts.sum()

    # Ensure lengths match for comparison
    min_length = min(len(expected_cdf), len(empirical_cdf))
    expected_counts = expected_counts[:min_length]
    expected_cdf = expected_cdf[:min_length]
    empirical_counts = empirical_counts[:min_length]
    empirical_cdf = empirical_cdf[:min_length]

    # Step 5: Perform Kolmogorov-Smirnov test
    ks_stat, p_value = kstest(empirical_cdf, expected_cdf)

    # Calculate the estimated percentage of cells expected to have at least one lentivirus (1 - P(X=0))
    percent_with_lentivirus = (1 - np.exp(-lambda_estimate)) * 100

    # Compile results into a dictionary
    results = {
        'Observed Mean': observed_mean,
        'Estimated Lambda': lambda_estimate,
        'K-S Statistic': ks_stat,
        'P-Value': p_value,
        'Percent with Lentivirus': percent_with_lentivirus
    }

    return empirical_counts, expected_counts, empirical_cdf, expected_cdf, results

def generate_clone_barcodes(filtered_df, max_lineage_count):
    """
    Generates clone barcodes and assigns clone IDs to each cell based on the filtered lineage barcodes.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        DataFrame containing the filtered lineage barcode data with columns 'cell_barcode' and 'corrected_lineage_barcode'.
    max_lineage_count : int, optional
        Maximum lineage count to consider for determining clone barcodes. Cells with lineage counts above this value will be excluded from clone barcode analysis.
        If not specified, all lineage counts will be considered.

    Returns
    -------
    tuple
        A tuple containing:
        - cloneid_clonebarcode_df (pd.DataFrame): DataFrame with columns 'cloneID' and 'clone_barcode'.
        - cloneid_cellbarcode_df (pd.DataFrame): DataFrame with columns 'cloneID' and 'cell_barcode'.
    """

    column_names = ['cell_barcode', 'corrected_lineage_barcode']

    # Filter the DataFrame based on max_lineage_count if provided
    if max_lineage_count:
        clone_barcode_df = filtered_df[filtered_df['lineage_count'] < max_lineage_count]
        clone_barcode_df = clone_barcode_df[column_names]
    else:
        clone_barcode_df = filtered_df[column_names]

    # Determine clone barcode for each cell_barcode
    def determine_clone_barcode(group):
        if len(group) == 1:
            return group.iloc[0]
        else:
            return "_".join(sorted(group))

    clone_barcode_df['clone_barcode'] = clone_barcode_df.groupby('cell_barcode')['corrected_lineage_barcode'].transform(determine_clone_barcode)

    # Drop duplicates to have one row per cell_barcode with its clone_barcode
    unique_barcodes = clone_barcode_df[['cell_barcode', 'clone_barcode']].drop_duplicates()

    # Determine the number of cell_barcodes associated with each clone_barcode
    clone_counts = unique_barcodes['clone_barcode'].value_counts().reset_index()
    clone_counts.columns = ['clone_barcode', 'count']

    # Sort the clone_barcodes by the number of associated cell_barcodes and assign cloneID
    clone_counts = clone_counts.sort_values(by='count', ascending=False).reset_index(drop=True)
    clone_counts['cloneID'] = clone_counts.index + 1

    # Generate the two required DataFrames
    cloneid_clonebarcode_df = clone_counts[['cloneID', 'clone_barcode']]
    cloneid_cellbarcode_df = unique_barcodes.merge(cloneid_clonebarcode_df, on='clone_barcode')
    cloneid_cellbarcode_df = cloneid_cellbarcode_df[['cloneID', 'cell_barcode']].sort_values(by="cloneID")

    return cloneid_clonebarcode_df, cloneid_cellbarcode_df, clone_counts


def save_results(filebase, z_threshold, max_lineage_count,  df, filtered_df, cloneid_clonebarcode_df, cloneid_cellbarcode_df, qc_results):
    """
    Saves the results of the analysis to the specified output directory.

    Parameters
    ----------
    filebase : str
        Directory path and prefix to save the output files and plots. 
    df : pd.DataFrame
        Original DataFrame with added normalization, log transformation, and Z-score columns.
    filtered_df : pd.DataFrame
        DataFrame with noise filtered out based on the Z-score threshold.
    cloneid_clonebarcode_df : pd.DataFrame
        DataFrame with columns 'cloneID' and 'clone_barcode'.
    cloneid_cellbarcode_df : pd.DataFrame
        DataFrame with columns 'cell_barcode' and 'cloneID'.
    qc_results : dict
        Dictionary containing QC results including observed mean, estimated lambda, K-S statistic, p-value, and percentage of cells with at least one lentivirus.

    Returns
    -------
    None
    """
    if max_lineage_count:
        zcfilebase = f'{filebase}.z{z_threshold}_mlc{max_lineage_count}'
        zfilebase = f'{filebase}.z{z_threshold}'
    else:
        zfilebase = f'{filebase}.z{z_threshold}'
        zcfilebase = f'{filebase}.z{z_threshold}'

    # Create the output directory if it does not exist
    # os.makedirs(output_dir, exist_ok=True)

    # Save the DataFrames to TSV files
    df.to_csv(f'{filebase}.rewind_statistics.tsv', sep='\t', index=False)
    filtered_df.to_csv(f'{zfilebase}.rewind_filtered.tsv', sep='\t', index=False)
    cloneid_clonebarcode_df.to_csv(f'{zcfilebase}.cloneID_cloneBarcode.tsv', sep='\t', index=False)
    cloneid_cellbarcode_df.to_csv(f'{zcfilebase}.cloneID_cellBarcode.tsv', sep='\t', index=False)

    # Save the QC results to a text file
    with open(f'{zfilebase}.truncated_Poisson_fit_results.txt', 'w') as f:
        for key, value in qc_results.items():
            f.write(f"{key}: {value}\n")

def plot_lineage_barcode_qc(filebase, empirical_counts, expected_counts, expected_cdf, empirical_cdf):
    """
    Generates QC plots for lineage barcode analysis, comparing empirical and expected counts, fractions, and CDFs.

    Parameters
    ----------
    filebase : str
        Directory and prefix to save the QC plot.
    empirical_counts : np.ndarray
        Array of empirical counts of lineage barcodes per cell.
    expected_counts : np.ndarray
        Array of expected counts from the truncated Poisson distribution.
    expected_cdf : np.ndarray
        Cumulative distribution function (CDF) of the expected counts.
    empirical_cdf : np.ndarray
        CDF of the empirical counts.

    Returns
    -------
    None
    """

    # Create x-axis values
    x_values = np.arange(1, len(expected_counts) + 1)

    # Define plotting data
    plot_data = {
        1: {
            'y_expected': expected_counts,
            'y_empirical': empirical_counts,
            'ylabel': 'Counts',
            'title': 'Expected vs Empirical Counts'
        },
        2: {
            'y_expected': expected_counts / expected_counts.sum(),
            'y_empirical': empirical_counts / empirical_counts.sum(),
            'ylabel': 'Fraction',
            'title': 'Expected vs Empirical Fractions'
        },
        3: {
            'y_expected': expected_cdf,
            'y_empirical': empirical_cdf,
            'ylabel': 'CDF',
            'title': 'Expected vs Empirical CDFs'
        }
    }

    # Create a figure with 3 subplots
    plt.figure(figsize=(18, 6))

    for i in range(1, 4):
        plt.subplot(1, 3, i)
        plt.plot(x_values, plot_data[i]['y_expected'], label='Expected')
        plt.plot(x_values, plot_data[i]['y_empirical'], label='Empirical')
        plt.xlabel('Number of Lineages per Cell')
        plt.ylabel(plot_data[i]['ylabel'])
        plt.legend()
        plt.title(plot_data[i]['title'])

    plt.tight_layout()
    plt.savefig(f'{filebase}.truncated_Poisson_fit_plots.png')
    plt.show()

def plot_clone_barcode_qc(filebase, clone_counts):
    """
    Generates QC plots for clone barcode analysis, showing the distribution of the number of cells per clone.

    Parameters
    ----------
    filebase : str
        Directory and prefix to save the QC plot.
    clone_counts : pd.Series
        Series containing the counts of cells per clone.

    Returns
    -------
    None
    """

    data = clone_counts['count']
    plt.figure(figsize = (6,6))
    # Create the histogram plot with fractions (density)
    sns.histplot(data, bins=range(1, max(data) + 2), discrete=True, stat='density')

    # Get the current axis
    ax = plt.gca()

    # Set the x-ticks to go from 1 to the max value of the data, incremented by 1
    ax.set_xticks(range(1, max(data) + 1))

    # Optional: Set the x-tick labels to be the same as the x-ticks
    ax.set_xticklabels(range(1, max(data) + 1))

    # Set the x-axis label
    ax.set_xlabel("Number of Cells per Clone")

    # Set the y-axis label to reflect fractions
    ax.set_ylabel("Fraction of Clones")

    # Save the plot
    plt.tight_layout()
    plt.savefig(f'{filebase}.clone_barcode_distribution_qc_plot.png')
    plt.show()

def main():
    """
    Main function to execute the Lineage and Clone Barcode Analysis.

    This script performs two main tasks:
    1. Filters out noise in lineage barcode data by normalizing, log-transforming, and calculating Z-scores to identify significant deviations.
    2. Determines clone barcodes and assigns clone IDs to each cell based on the filtered lineage barcodes, and generates QC plots for both tasks.

    Steps:
    - Load input data.
    - Execute Task 1: Filter out noise.
    - Execute Task 1 QC: Test if the distribution fits a truncated Poisson distribution.
    - Execute Task 2: Determine clone barcodes and assign clone IDs.
    - Execute Task 2 QC: Plot the distribution of the number of cells per clone.
    - Save results to the specified output directory.
    - Plot and save QC plots for both tasks.
    """
    args = parse_args()

    # use as the directory and prefix of outputs
    filebase = args.input_file.rsplit(".", 1)[0]
    if args.max_lineage_count:
        zcfilebase = f'{filebase}.z{args.z_threshold}_mlc{args.max_lineage_count}'
        zfilebase = f'{filebase}.z{args.z_threshold}'
    else:
        zfilebase = f'{filebase}.z{args.z_threshold}'
        zcfilebase = f'{filebase}.z{args.z_threshold}'

    print(filebase)
    print(zfilebase)
    print(zcfilebase)


    # Load input data
    df = pd.read_csv(args.input_file, sep="\t")

    # Task 1: Filter out noise
    df, filtered_df = filter_noise(df, args.top_n, args.z_threshold, args.min_group_size)
    # Task 1 QC: Test if the distribution of total lineage barcodes per cell fits a truncated Poisson distribution
    lineage_count_df = filtered_df[['cell_barcode', 'lineage_count']].drop_duplicates().reset_index(drop=True)
    empirical_counts, expected_counts, empirical_cdf, expected_cdf, results = test_truncated_poisson_distribution(lineage_count_df, 'lineage_count')
    # Plot QC plots for Task 1
    plot_lineage_barcode_qc(zfilebase, empirical_counts, expected_counts, expected_cdf, empirical_cdf)


    # Task 2: Determine clone barcodes and assign clone IDs
    cloneid_clonebarcode_df, cloneid_cellbarcode_df, clone_counts = generate_clone_barcodes(filtered_df, args.max_lineage_count)
    # Plot QC plots for Task 2
    # clone_counts = cloneid_clonebarcode_df['clone_barcode'].value_counts()
    plot_clone_barcode_qc(zcfilebase, clone_counts)


    # Save results
    save_results(filebase, args.z_threshold, args.max_lineage_count, df, filtered_df, cloneid_clonebarcode_df, cloneid_cellbarcode_df, results)


if __name__ == "__main__":
    main()
