# %%
# Declare the important variables with their contents
bases = ["A", "T", "C", "G"]
gc_bases = ["C", "G"]
complement_bases = {
    "A": "T", "T": "A", "G": "C", "C": "G"
}
motifs = ["ATG", "TAA", "TAG", "TGA"]


# %%
# Check if the input sequence is valid
def is_valid(sequence):
    for base in sequence:
        if base not in bases:
            return False
    
    # If all the bases are in ATCG
    return True


# %%
# Print the length of the input sequence
def sequence_properties(sequence):
    sequence_length = len(sequence)
    print(f"Length: {sequence_length} base pairs.\n")


# %%
def count_bases(sequence):
    # Create a dictionary to be filled with each base and its occurrence percentage
    count = {}
    # Calculate each base percent and append it to the dictionary
    for base in sequence:
        base_percent = round(sequence.count(base)/len(sequence) *100, 2)
        count[base] = f"{base_percent}%"
    # Print the components of the dictionary
    for base, percent in count.items():
        print(f"{base}: {percent}")


# %%
def gc_content(sequence):
    # Calculate the sequence's gc content
    gc_count = 0
    for base in sequence:
        if base in gc_bases:
            gc_count += 1
    gc_percent = round((gc_count/len(sequence)) * 100, 2)
    return gc_percent


# %%
def reverse_complement(sequence):
    # Create an empty string for the reverse complement
    complement = ""

    # For each base in the sequence, return its cmplement using the earlier defined complement dictionary
    for base in sequence:
        complement += complement_bases[base]

    # Reverse the complement string and print it
    rev_complement = complement[::-1]
    print(f"Reverse Complement: {rev_complement}")


# %%
def find_motif(sequence):

    # Looping through every three position
    for i in range(0, len(sequence), 3):

        # Checks if a base and the next two (a codon) is a motif, and if it's a start or stop codon
        if sequence[i:i+3] in motifs:
            if sequence[i:i+3] == "ATG":
                print(f"Start codon found at position {i+1}")
            else:
                print(f"Stop codon found at position {i+1}")

# %%
#Create a function that runs all the other functions if the sequence is valid
def analyze_sequence():
    # Take in a sequence from the user
    sequence = input("Please enter a valid sequence:")
    sequence = sequence.upper()

    print("DNA Sequence Analysis\n======================================================")
    print(f"DNA Sequence:{sequence}\n")

    # Runs all the other funtions if the sequence is valid
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

# Run the aggregator function
analyze_sequence()

