# %%
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
from Bio.Data import CodonTable
from itertools import product
from collections import Counter
from datetime import datetime
import os
import glob
import sys

start_codon = "ATG"
stop_codon = ["TAA", "TAG", "TGA"]
input_files = glob.glob("test_data/test*.fasta*")
global_codon_counts = Counter()
total_num_of_codons = 0
each_codon_freq = {}
results = []
now = datetime.now()
formatted = now.strftime("%Y-%m-%d %H:%M:%S")

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

def calculate_frequencies(global_codon_counts, total_num_of_codons):
    total_num_of_codons = sum(global_codon_counts.values())
    for codon in global_codon_counts.keys():
        codon_freq = round(global_codon_counts[codon]/total_num_of_codons, 3)
        each_codon_freq[codon] = codon_freq

    return each_codon_freq 

def save_csv(global_codon_counts, frequencies):
    with open("codon_frequency.csv", "w") as file:
        file.write("Codon,Count,Frequency,Percentage\n")
        for codon, value in global_codon_counts.items():
            file.write(f"{codon},{value},{frequencies[codon]},{frequencies[codon] * 100}%\n")

def generate_report():
    with open("summary.txt", "w") as file:
        file.write("Codon Usage Analysis Report\n")
        file.write(f"{'=' * 30}\n")
        file.write(f"Generated on: {formatted}\n")
        file.write(f"Input files: \n")
        file.write(f"   \n")
        

for file in input_files:
    print(f"This is for the {file.upper()} file")
    file_codon_counts = Counter()
    records = SeqIO.parse(file, "fasta")
    for record in records:
        sequence = record.seq
        seq_id = record.id
        is_valid, message = validate_cds(sequence, seq_id)  
        print(f"{seq_id}: {message}")

        if is_valid:
            codons = extract_codons(sequence)
            codon_counts = count_codons(codons)
            file_codon_counts += codon_counts
            global_codon_counts += codon_counts

            results.append({
                "id": seq_id,
                "length": len(sequence),
                "gc_content": round(gc_fraction(sequence) * 100, 2),
                "codon_counts": codon_counts,            
            })

    print(file_codon_counts)

print(f"The global codon counts: {global_codon_counts}")

start_count, stop_count, taa_count, tag_count, tga_count = identify_special_codons(global_codon_counts)
print(f"There are {start_count} start codons and {taa_count} TAA's, {tag_count} TAG's, and {tga_count} TGA's, making {stop_count} stop codons in total")

frequencies = calculate_frequencies(global_codon_counts, total_num_of_codons)
print(frequencies)

