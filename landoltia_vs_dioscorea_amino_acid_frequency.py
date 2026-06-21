!pip install biopython
!pip install matplotlib

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import matplotlib.pyplot as plt
from collections import Counter

Entrez.email = "sahana.arvind.2010@gmail.com"

# load nucleotide sequence from Landoltia punctata
handle = Entrez.efetch(db="nucleotide", id="JZ990160.1", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

seq10 = record.seq

# translate nucleotide to protein
translated_protein = record.seq.translate(to_stop=True)
print(f"Translated protein (default): {translated_protein}")

# load protein sequence from Dioscorea alata
handle = Entrez.efetch(db="protein", id="KAH7657379.1", rettype="fasta", retmode="text")
prot_record = SeqIO.read(handle, "fasta")
handle.close()

protein_seq = prot_record.seq
print(f"Protein sequence: {protein_seq}")

# count amino acids, standard amino acids only
def count_amino_acids(sequence):
    sequence = str(sequence).upper()
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    filtered_seq = [aa for aa in sequence if aa in valid_aas]
    return Counter(filtered_seq)

counts_translated = count_amino_acids(translated_protein)
counts_protein = count_amino_acids(protein_seq)

# Calculate total valid amino acids for normalization
total_translated = sum(counts_translated.values())
total_protein = sum(counts_protein.values())

# list standard amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

# --- MODIFIED SECTION: Calculate percentages instead of raw counts ---
# (count / total_valid_aas) * 100
pct_translated_list = [(counts_translated.get(aa, 0) / total_translated) * 100 for aa in amino_acids]
pct_protein_list = [(counts_protein.get(aa, 0) / total_protein) * 100 for aa in amino_acids]

# Landoltia punctata graph
plt.figure(figsize=(10, 5))
plt.bar(amino_acids, pct_translated_list, color="darkblue")

plt.xlabel("Amino Acid")
plt.ylabel("Percentage (%)") # Updated label
plt.title("Amino Acid Composition (%): 19SA1.24 (t-snare protein-like mRNA) Protein Translation")

plt.tight_layout()
plt.savefig("landoltia_amino_acid_percentage.png", dpi=300)
plt.show()

# Dioscorea alata graph
plt.figure(figsize=(10, 5))
plt.bar(amino_acids, pct_protein_list, color="darkblue")

plt.xlabel("Amino Acid")
plt.ylabel("Percentage (%)") # Updated label
# FIXED EDITOR COMMENT 3: Ensured title strictly matches the fetched accession KAH7657379.1
plt.title("Amino Acid Composition (%): KAH7657379.1 t-snare protein [Dioscorea alata]")

plt.tight_layout()
plt.savefig("dioscorea_amino_acid_percentage.png", dpi=300)
plt.show()

# print the percentages for your verification
print("Translated protein percentages:")
for aa, pct in zip(amino_acids, pct_translated_list):
    if pct > 0: print(f"{aa}: {pct:.2f}%")

print("\nProtein sequence percentages:")
for aa, pct in zip(amino_acids, pct_protein_list):
    if pct > 0: print(f"{aa}: {pct:.2f}%")
