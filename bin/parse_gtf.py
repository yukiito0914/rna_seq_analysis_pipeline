#!/usr/bin/env python

# You can refer to the help manual by `python genome_stats.py -h`

import argparse
import re

# Initialize argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help='a GTF file containing genomic annotations', dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file where we will write gene_id and gene_name mappings (tab-delimited)', dest="output", required=True)
args = parser.parse_args()

# Initialize a set to store unique (gene_id, gene_name) pairs
seen_pairs = set()

# Open input and output files
with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    # Write header line
    outfile.write("gene_id\tgene_name\n")
    
    # Process each line in the input file
    for line in infile:
        line = line.strip()
        # Skip empty lines or comments
        if not line or line.startswith("#"):
            continue
        
        # Split line by tabs; a valid GTF line should have at least 9 fields
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        
        # The 9th field contains attributes (e.g., gene_id "ENSG..."; gene_name "XYZ"; ...)
        attributes = fields[8]
        
        # Extract gene_id and gene_name using regular expressions
        m_gene_id = re.search(r'gene_id "([^"]+)"', attributes)
        m_gene_name = re.search(r'gene_name "([^"]+)"', attributes)
        
        # If both gene_id and gene_name were found, process them
        if m_gene_id and m_gene_name:
            gene_id = m_gene_id.group(1)
            gene_name = m_gene_name.group(1)
            pair = (gene_id, gene_name)
            # Only write if this pair hasn't been seen before
            if pair not in seen_pairs:
                seen_pairs.add(pair)
                outfile.write(f"{gene_id}\t{gene_name}\n")
