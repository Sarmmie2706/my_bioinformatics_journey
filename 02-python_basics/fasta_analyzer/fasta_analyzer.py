# %%
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight

# %%
records = list(SeqIO.parse("test_genes.fasta", "fasta"))
record_list = []
for record in records:
    record_dict = {}
    record_dict["id"] = record.id
    record_dict["length"] = len(record.seq)
    record_dict["gc_fraction"] = round(gc_fraction(record.seq)*100, 2)
    record_dict["start_count"] = record.seq.count("ATG")
    record_dict["stop_count"] = record.seq.count("TAG") +record.seq.count("TAA") 
    + record.seq.count("TGA")
    record_dict["molecular_weight"] = round(molecular_weight(record.seq, "DNA"), 2)
    record_list.append(record_dict)
print(record_list)


# %%
print(f"{'ID':<10}{'Length':<10}{'GC%':<10}{'Start':<8}{'Stop':<8}{'Mol.Weight':<12}")
print("-" * 56)
for record in record_list:
    print(f"{record['id']:<10}{record['length']:<10}{record['gc_fraction']:<10}{record['start_count']:<8}{record['stop_count']:<8}{record['molecular_weight']:<12}")


# %%
print("Summary Statistics")
print("=" * 65)
print(f"Total Sequence: {len(record_list)}")
total_length = 0
for record in record_list:
    total_length += record["length"]
ave_length = round((total_length/len(record_list)),2)
print(f"Average length: {ave_length}bp")

max_length = 0
max_id = ""
for record in record_list:
    if record["length"] > max_length:
        max_length = record["length"]
        max_id = record["id"]
print(f"Longest: {max_length} bp ({max_id})")

min_length = 160000000000
min_id = ""
for record in record_list:
    if record["length"] < min_length:
        min_length = record["length"]
        min_id = record["id"]
print(f"Longest: {min_length} bp ({min_id})")

total_gc = 0
for record in record_list:
    total_gc += record["gc_fraction"]
ave_gc = round((total_gc/len(record_list)),2)
print(f"Average GC%: {ave_gc}bp")


# %%
with open("analyze_results.csv", "w") as file:
    file.write(f"ID','Length','GC%','Start','Stop','Mol.Weight'\n")
    for record in record_list:
        file.write(f"{record['id']},{record['length']},{record['gc_fraction']},{record['start_count']},{record['stop_count']},{record['molecular_weight']}\n")

# %%
with open("summary.txt", "w") as file:
    file.write("Summary Statistics\n")
    file.write(f"{'=' * 65}\n")
    file.write(f"Total Sequence: {len(record_list)}\n")
    file.write(f"Average length: {ave_length}bp\n")
    file.write(f"Average GC%: {ave_gc}bp\n")
    file.write(f"Longest: {max_length} bp ({max_id})\n")
    file.write(f"Shortest: {min_length} bp ({min_id})\n")



