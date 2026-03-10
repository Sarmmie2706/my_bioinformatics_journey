# %%
sequence = "gtatagatgcccgctaacgtctaggcagtacggattccagtacac"
bases = ["A", "T", "C", "G"]
lower_bases = [base.lower() for base in bases]
gc_bases = ["C", "G", "c", "g"]
complement_bases = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "a": "t", "t": "a", "g": "c", "c": "g"
}
motifs = ["ATG", "TAA", "TAG", "TGA"]
lower_motifs = [motif.lower() for motif in motifs]

# %%
def is_valid(sequence):
    for base in sequence:
        if base not in bases and base not in lower_bases:
            print("This is not a valid DNA sequence.")
            return False
    
    # If all the bases are in ATCG
    return True
is_valid(sequence)  

# %%
def sequence_properties(sequence):
    sequence_length = len(sequence)
    print(f"The sequence is {sequence_length} base pairs long\n")


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
gc_content(sequence)

# %%
def reverse_complement(sequence):
    complement = ""
    for base in sequence:
        complement += complement_bases[base]
    rev_complement = complement[::-1]
    return rev_complement
reverse_complement(sequence)

# %%
def find_motif(sequence):
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] in motifs or sequence[i:i+3] in lower_motifs:
            return i+1
find_motif(sequence)

# %%
def analyze_sequence():
    # sequence = input("Please enter a valid sequence:")
    sequence = "gtatagatgcccgctaacgtctaggcagtacggattccagtacac"

    print("DNA Sequence Analysis\n======================")
    print(f"DNA Sequence:{sequence}\n")
    if is_valid(sequence) == True:
        print("This is a valid DNA sequence.\n")
        sequence_properties(sequence)
        
        print("Base Counts")
        count_bases(sequence)

        print("\n")
        print(f"It's GC content is {gc_content(sequence)}%")




analyze_sequence()

# %%


# %%
