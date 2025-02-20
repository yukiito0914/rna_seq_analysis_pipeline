#!/usr/bin/env python

# You can refer to the help manual by `python genome_stats.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help='a GTF file containing a genomic annotations', dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file where we will gene_id and gene_name mappings (tab-delimited)', dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()


import re

# Open the input GTF file and the output file
with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    # Write the header line for the output file
    outfile.write("gene_id\tgene_name\n")
    
    # Process each line in the input file
    for line in infile:
        line = line.strip()
        # Skip comment lines that start with '#' or empty lines
        if not line or line.startswith("#"):
            continue
        
        # Split the line by tabs; a valid GTF line should have at least 9 fields
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        
        # The 9th field contains attributes (e.g., gene_id "ENSG..."; gene_name "XYZ"; ...)
        attributes = fields[8]
        
        # Use a regular expression to extract the gene_id from the attributes field
        m_gene_id = re.search(r'gene_id "([^"]+)"', attributes)
        # Use a regular expression to extract the gene_name from the attributes field
        m_gene_name = re.search(r'gene_name "([^"]+)"', attributes)
        
        # If both gene_id and gene_name were found, write them to the output file
        if m_gene_id and m_gene_name:
            gene_id = m_gene_id.group(1)
            gene_name = m_gene_name.group(1)
            outfile.write(f"{gene_id}\t{gene_name}\n")