#!/usr/bin/env python

import argparse
import os
import pandas as pd

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        description="Combine count outputs from each sample into a single counts matrix using pandas"
    )
    # Accept a list of input count file paths; this allows bracket notation without quotes
    parser.add_argument("-i", "--input", nargs='+', help="List of input count file paths. Example: [A.txt, B.txt, C.txt]", required=True)
    # Output file path for the combined counts matrix
    parser.add_argument("-o", "--output", help="Output file for the combined counts matrix", required=True)
    args = parser.parse_args()

    files = args.input

    # Check if the input is provided in bracket notation split across arguments
    if len(files) > 1 and files[0].startswith('[') and files[-1].endswith(']'):
        # Join the tokens into a single string, then remove the surrounding brackets
        input_str = " ".join(files)
        input_str = input_str.strip()[1:-1]
        # Split by comma and remove extra whitespace
        files = [f.strip() for f in input_str.split(',') if f.strip()]

    if not files:
        print("No input files provided.")
        exit(1)

    # List to store DataFrames for each sample
    dfs = []
    # List to store sample names (derived from file names)
    sample_names = []

    # Process each file in the list
    for filepath in files:
        # Extract sample name from file name (removing the extension)
        sample_name = os.path.basename(filepath).split('.')[0]
        sample_names.append(sample_name)

        # Read the count file into a DataFrame (tab-delimited with headers, e.g., 'gene' and 'count')
        df = pd.read_csv(filepath, sep='\t')
        # Rename the 'count' column to the sample name
        df = df.rename(columns={"count": sample_name})
        # Set the 'gene' column as the index for merging on gene names
        df = df.set_index("gene")
        dfs.append(df)

    # Merge all DataFrames on the gene index using an outer join to include all genes
    combined_df = pd.concat(dfs, axis=1, join="outer")
    # Fill missing values (i.e., genes not found in a sample) with 0
    combined_df = combined_df.fillna(0)
    # Reset the index so that gene names become a regular column
    combined_df = combined_df.reset_index()

    # Sort the sample columns alphabetically for consistent order
    sorted_samples = sorted(sample_names)
    combined_df = combined_df[["gene"] + sorted_samples]

    # Write the combined counts matrix to the output file as a tab-delimited text file
    combined_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()