# %%
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
from Bio.Data import CodonTable
from itertools import product
from collections import Counter
import os
import glob
import sys

start_codon = "ATG"
stop_codon = ["TAA", "TAG", "TGA"]
input_files = glob.glob("test_data/test*.fasta*")

def validate_cds(sequence, seq_id):
    if str(sequence[:3]) == start_codon and str(sequence[-3:]) in stop_codon and len(sequence)%3 == 0:
        return True, "Valid CDS"
    else:
        return False, "This is not a valid coding sequence"

def extract_codons(sequence):
    sequence_list = []
    for i in range(0, len(sequence), 3):
        sequence_list.append(str(sequence[i:i+3]))
    return sequence_list

def count_codons(codons):
    return Counter(codons)

def identify_special_codons(global_codon_counts):
    start_count = 0
    stop_count = 0
    taa_count = 0
    tga_count = 0
    tag_count = 0
    for key in global_codon_counts.keys():
        if key == "ATG":
            start_count += global_codon_counts[key]
        elif key == "TAA":
            taa_count += global_codon_counts[key]
            stop_count += global_codon_counts[key]
        elif key == "TGA":
            tga_count += global_codon_counts[key]
            stop_count += global_codon_counts[key]
        elif key == "TAG":
            tag_count += global_codon_counts[key]
            stop_count += global_codon_counts[key]
    return start_count, stop_count, taa_count, tag_count, tga_count

global_codon_counts = Counter()

for file in input_files:
    print(f"This is for the {file.upper()} file")
    file_codon_counts = Counter()
    records = SeqIO.parse(file, "fasta")
    for record in records:
        sequence = record.seq
        seq_id = record.id
        is_valid, message = validate_cds(sequence, seq_id)  # call here
        print(f"{seq_id}: {message}")

        if is_valid:
            codons = extract_codons(sequence)
            codon_counts = count_codons(codons)
            file_codon_counts += codon_counts
            global_codon_counts += codon_counts
    print(file_codon_counts)

print(f"The global codon counts: {global_codon_counts}")

start_count, stop_count, taa_count, tag_count, tga_count = identify_special_codons(global_codon_counts)
print(f"There are {start_count} start codons and {taa_count} TAA's, {tag_count} TAG's, and {tga_count} TGA's, making {stop_count} stop codons in total")



