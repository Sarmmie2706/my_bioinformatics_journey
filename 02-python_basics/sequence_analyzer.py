# %%
bases = ["A", "T", "C", "G"]
gc_bases = ["C", "G"]
complement_bases = {
    "A": "T", "T": "A", "G": "C", "C": "G"
}
motifs = ["ATG", "TAA", "TAG", "TGA"]


# %%
def is_valid(sequence):
    for base in sequence:
        if base not in bases:
            return False
    
    # If all the bases are in ATCG
    return True


# %%
def sequence_properties(sequence):
    sequence_length = len(sequence)
    print(f"Length: {sequence_length} base pairs.\n")


# %%
def count_bases(sequence):
    count = {}
    for base in sequence:
        base_percent = round(sequence.count(base)/len(sequence) *100, 2)
        count[base] = f"{base_percent}%"
    for base, percent in count.items():
        print(f"{base}: {percent}")


# %%
def gc_content(sequence):
    gc_count = 0
    for base in sequence:
        if base in gc_bases:
            gc_count += 1
    gc_percent = round((gc_count/len(sequence)) * 100, 2)
    return gc_percent


# %%
def reverse_complement(sequence):
    complement = ""
    for base in sequence:
        complement += complement_bases[base]
    rev_complement = complement[::-1]
    print(f"Reverse Complement: {rev_complement}")


# %%
def find_motif(sequence):
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] in motifs:
            if sequence[i:i+3] == "ATG":
                print(f"Start codon found at position {i+1}")
            else:
                print(f"Stop codon found at position {i+1}")

# %%
def analyze_sequence():
    sequence = input("Please enter a valid sequence:")
    sequence = sequence.upper()

    print("DNA Sequence Analysis\n======================================================")
    print(f"DNA Sequence:{sequence}\n")
    if is_valid(sequence) == True:
        print("Valid DNA sequence.\n")
        sequence_properties(sequence)
        
        print("Base Counts")
        count_bases(sequence)

        print("\n")
        print(f"GC content: {gc_content(sequence)}%")

        print("\n")
        reverse_complement(sequence)

        print("\n")
        find_motif(sequence)

        print("\n")
        print("Analysis Complete!")

    else:
        print("Please enter a valid DNA sequence")

analyze_sequence()


# %%
