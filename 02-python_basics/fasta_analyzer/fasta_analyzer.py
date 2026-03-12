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
