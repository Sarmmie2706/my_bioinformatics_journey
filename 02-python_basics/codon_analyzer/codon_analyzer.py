# %%
# Import the required modules
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

# Set all the needed global variables that will be used later inside functions
start_codon = "ATG"
stop_codon = ["TAA", "TAG", "TGA"]

# This searches the filesystem and adds the files as a list
input_files = glob.glob("test_data/test*.fasta*")

global_codon_counts = Counter()
total_num_of_codons = 0
each_codon_freq = {}
results = []
now = datetime.now()
formatted = now.strftime("%Y-%m-%d %H:%M:%S")
total_sequences = 0
valid_sequences = 0

# Confirm if a sequence is a valid coding sequence
def validate_cds(sequence, seq_id):
    # Checks if the first three and last three bases are start and stop codons
    # And also if the codon length is divisible by 3
    if str(sequence[:3]) == start_codon and str(sequence[-3:]) in stop_codon and len(sequence)%3 == 0:
        return True, "Valid CDS"
    else:
        return False, "This is not a valid coding sequence"

# Extracts all the codons and compiles in a list
def extract_codons(sequence):
    sequence_list = []
    for i in range(0, len(sequence), 3):
        sequence_list.append(str(sequence[i:i+3]))
    return sequence_list

# Counts all the codons
def count_codons(codons):
    return Counter(codons)

# Takes in the counter to count the number of start codons, 
# total and individual stop codons and returns the figures
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

# Calculating each codon's frequency
def calculate_frequencies(global_codon_counts, total_num_of_codons):
    # Add up all values from the global codon counts' values
    total_num_of_codons = sum(global_codon_counts.values())

    # Divide each codon's individual count by the total to get its frequency (percent)
    # Each frequency is appended to a dictionary
    for codon in global_codon_counts.keys():
        codon_freq = round(global_codon_counts[codon]/total_num_of_codons, 3)
        each_codon_freq[codon] = codon_freq

    # The function returns a dictionary with the codon and its frequency
    return each_codon_freq 

# Save important results to a csv file
def save_csv(global_codon_counts, frequencies):
    with open("codon_frequency.csv", "w") as file:
        file.write("Codon,Count,Frequency,Percentage\n")
        for codon, value in global_codon_counts.items():
            file.write(f"{codon},{value},{frequencies[codon]},{round(frequencies[codon] * 100, 3)}%\n")

# Generating a summary report
def generate_report():
    with open("summary.txt", "w") as file:
        file.write("Codon Usage Analysis Report\n")
        file.write(f"{'=' * 30}\n")
        file.write(f"Generated on: {formatted}\n")
        for f in input_files:
            file.write(f"   - {f}\n")
        file.write(f"   \n")
        file.write("Summary\n")
        file.write(f"{'-' * 8}\n")
        file.write(f"Total sequences analyzed: {total_sequences}\n")
        file.write(f"Valid coding sequences: {valid_sequences}\n")
        file.write(f"{'-' * 8}\n")
        file.write(f"Start codons\n")
        file.write(f"{start_codon}: {start_count}\n")
        file.write(f"{'-' * 8}\n")
        file.write(f"Stop codons\n")
        file.write(f"{stop_codon[0]}: {taa_count}\n")
        file.write(f"{stop_codon[1]}: {tag_count}\n")
        file.write(f"{stop_codon[2]}: {tga_count}\n")
        file.write(f"{'-' * 8}\n")
        file.write("Most common codons\n")

        # most_common and least_common return a list of tuple containing (codon, count)
        # that can be looped through
        for a,b in most_common:
            file.write(f"{a}: {b}\n")
        file.write(f"{'-' * 8}\n")
        file.write("Least common codons\n")
        for a,b in least_common:
            file.write(f"{a}: {b}\n")
        file.write(f"{'-' * 8}\n")
        file.write("Files created\n")
        file.write("- codon_frequency.csv\n")
        file.write("- analysis_report.txt\n")

# Analysis for each file
for file in input_files:
    print(f"Analysing {file.upper()}")
    file_codon_counts = Counter()
    records = SeqIO.parse(file, "fasta")
    for record in records:
        total_sequences += 1
        sequence = record.seq
        seq_id = record.id
        is_valid, message = validate_cds(sequence, seq_id)  
        print(f"{seq_id}: {message}")

        # If a sequence is valid, calculate its codon counts from the extract_codons and count_codons
        # function and increment the file and global codon counts based on the results
        if is_valid:
            valid_sequences += 1
            codons = extract_codons(sequence)
            codon_counts = count_codons(codons)
            file_codon_counts += codon_counts
            global_codon_counts += codon_counts

            # Each sequence's results is also appended to the result variable
            results.append({
                "id": seq_id,
                "length": len(sequence),
                "gc_content": round(gc_fraction(sequence) * 100, 2),
                "codon_counts": codon_counts,            
            })
    # Print each file's codon counts
    print(file_codon_counts)

# Print the global codon counts
print(f"The global codon counts: {global_codon_counts}")

# The values returned from the identify_special_ codons functions are collected into variables
# These variables are eventually used in the generate_summary function
start_count, stop_count, taa_count, tag_count, tga_count = identify_special_codons(global_codon_counts)
print(f"There are {start_count} start codons and {taa_count} TAA's, {tag_count} TAG's, and {tga_count} TGA's, making {stop_count} stop codons in total")

# The each_codon_freq dictionary that is returned from the function is assigned to the frequencies variable
# This is later used in the save_csv function
frequencies = calculate_frequencies(global_codon_counts, total_num_of_codons)
print(frequencies)

# These are the codons with the highest and least count respectively
most_common = global_codon_counts.most_common(5)
least_common = global_codon_counts.most_common()[-5:]

# Calling the functions
save_csv(global_codon_counts, frequencies)
generate_report()




# %%
