#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Filter gene count data based on a threshold.")
parser.add_argument("-i", "--input", required=True, help="Input file (CSV, TSV, or TXT format)")
parser.add_argument("-t", "--threshold", type=float, required=True, help="Threshold for filtering")
parser.add_argument("-o", "--output", default="filtered_counts.csv", help="Output file name (default: filtered_counts.csv)")
parser.add_argument("-s", "--summary", default="summary.txt", help="Summary output file (default: summary.txt)")
args = parser.parse_args()

# Function to determine file separator
def guess_separator(file):
    with open(file, 'r') as f:
        first_line = f.readline()
    
    if ',' in first_line:
        return ','
    elif '\t' in first_line:
        return '\t'
    else:
        return '\t'  # Default to tab for .txt files

sep = guess_separator(args.input)

# Read data
df = pd.read_csv(args.input, sep=sep, header=0, index_col=0, dtype=str)

# Convert all columns to numeric
df = df.apply(pd.to_numeric, errors='coerce')

# Check if there are NA values introduced
if df.isna().any().any():
    print("Warning: NA values were introduced during conversion. Check input file format.")
    print(df.head())

# Remove rows with NA values
df = df.dropna()

# Filter genes based on the threshold
keep = df.sum(axis=1) > args.threshold
filtered_counts = df[keep]

# Save the filtered data to a file
filtered_counts.to_csv(args.output, sep=sep)

# Save summary to a text file
num_filtered = keep.sum()
num_removed = len(keep) - num_filtered
summary_text = f"Filtering threshold: {args.threshold}\nFiltered genes: {num_filtered}\nRemoved genes: {num_removed}\n"
with open(args.summary, "w") as f:
    f.write(summary_text)