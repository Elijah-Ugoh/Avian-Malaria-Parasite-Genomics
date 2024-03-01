# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:37:48 2024

@author: Elijah

Usage: exclude_bird_contigs.py scaffold.txt gffParse.faa > [provide_output_file_name]

Approach:
This script redas in the scaffold file containing the list of contigs deriving from bird genome
Then, it reads in the filtered fasta file obtained from the gene prediction
It runs a step-wise match by searching for each contig name in the fasta sequence headers. 
Next, it exclude the header and corresponding sequences that match contigs contained in the pasrsed contig list.
Finally, the scripts prints the result to standatrd output. This result is the gfiltered fasta file.   

"""
# Read the list of contigs to exclude
import sys

try:
    contiglist = sys.argv[1]
    fastafile = sys.argv[2]

except IndexError:
    print("Usage: exclude_bird_contigs.py scaffold.txt gffParse.faa")
    sys.exit()

exclude_contigs = set()  # create a set to hold the excluded contigs
with open(contiglist, "r") as exclude_file:
    for line in exclude_file:
        contig_id = line.strip()  # removing all white spaces from each contig name
        exclude_contigs.add(contig_id)  # add the contig to the set above

# Read the fasta file, excluding contigs present in the exclude list
current_contig = None
exclude_sequence = False
with open(fastafile, "r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            # Extract contig ID from the fasta header
            current_contig = line.split()[0][1:]  # Remove the leading '>'
            # Check if the contig should be excluded
            exclude_sequence = current_contig in exclude_contigs
            # Print header to the standard output if the sequence is not excluded
            if not exclude_sequence:
                print(line, end="")
        else:
            # Print sequence to the standard output if the sequence is not excluded
            if not exclude_sequence:
                print(line, end="")

