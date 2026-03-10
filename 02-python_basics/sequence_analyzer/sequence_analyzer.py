# %%
sequence = input("Please enter a valid DNA sequence:")
bases = ["A", "T", "C", "G"]
lower_bases = [base.lower() for base in bases]
gc_bases = ["C", "G", "c", "g"]
complement_bases = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "a": "t", "t": "a", "g": "c", "c": "g"
}
print(sequence)

# %%
def is_valid(sequence):
    for base in sequence:
        if base not in bases and base not in lower_bases:
            return "This is not a valid DNA sequence."
    
    # If all the bases are in ATCG
    return "This is a valid DNA sequence."

print(is_valid(sequence))    
    
# %%
def count_bases(sequence):
    count = {}
    for base in sequence:
        count[base] = sequence.count(base)
    return count
count_bases(sequence)

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
