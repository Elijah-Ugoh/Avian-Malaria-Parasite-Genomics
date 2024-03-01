# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 08:55:29 2024

@author: Elijah
Approach:
    This python script takes multiple protein fasta files and a ortholog (busco) table file a  as arguments
    It searches for each ortholog in each ortholog group in the fasta file where it is present, and then 
    saves the orthologs to a unique fasta file. 
    It produces as many fasta files as there are unique ortholog groups

Usage: python script.py busco_table.txt fasta_file1.txt fasta_file2.txt ...
"""
import csv
import sys

# Function to read the table of busco genes
def read_busco_table(busco_file):
    row_num = 0
    busco_genes = {}
    with open(busco_file, 'r') as csvfile:
        busco_groups = csv.reader(csvfile, delimiter='\t')
        next(busco_groups)  # Skip the header row
        for row in busco_groups:
            row_num += 1
            gene_ids = row[3:11]  # Assuming gene IDs are in columns 4 to 11
            busco_genes[row_num] = gene_ids
    return busco_genes

# Function to read a fasta file and return a dictionary with gene IDs as keys and sequences as values
def read_fasta_file(fasta_file):
    fasta_data = {}
    with open(fasta_file, 'r') as file:
        gene_id = None
        sequence = ""
        for line in file:
            if line.startswith('>'):
                if gene_id is not None:
                    fasta_data[gene_id] = sequence
                gene_id = line.strip()[1:]
                sequence = ""
            else:
                sequence += line.strip()
        if gene_id is not None:
            fasta_data[gene_id] = sequence
    return fasta_data

# Function to process busco data and fasta files
def process_busco_and_fasta(busco_genes):
    species = ['ht', 'pb', 'pc', 'pf', 'pk', 'pv', 'py', 'tg']
    gene_files = {species[i]: read_fasta_file(f'{species[i]}.faa') for i in range(len(species))}
    
    # Next, iterate through the busco groups and gene ids saved in the busco_genes dictionary
    for busco_group, gene_ids in busco_genes.items():
        with open(f'busco_{busco_group}.faa', "w") as busco_files:
            # Iterate through each busco group and write the orthologs to a unique file
            for species_index, gene in enumerate(gene_ids):
                # Get the species identifier from the list
                species_identifier = species[species_index]
                
                # Write the sequence ID with species identifier
                busco_files.write(f">{species_identifier}\n{gene_files[species_identifier][gene]}\n")
# Specify command line inputs            
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py busco_table.txt gene_file1.txt gene_file2.txt ...")
        sys.exit(1)
    # Get the busco table file and gene files from command-line arguments
    busco_genes = read_busco_table(sys.argv[1])
    read_fasta_files = sys.argv[2:]
    process_busco_and_fasta(busco_genes)