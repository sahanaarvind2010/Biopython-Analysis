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
handle = Entrez.efetch(db="protein", id="KAH7659748.1", rettype="fasta", retmode="text")
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

# list standard amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

counts_translated_list = [counts_translated.get(aa, 0) for aa in amino_acids]
counts_protein_list = [counts_protein.get(aa, 0) for aa in amino_acids]

# Landoltia punctata graph
plt.figure(figsize=(10, 5))
plt.bar(amino_acids, counts_translated_list, color="darkblue")

plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.title("Amino Acid Frequency: 19SA1.24 (t-snare protein-like mRNA) Protein Translation")

plt.tight_layout()
plt.savefig("landoltia_amino_acid_frequency.png", dpi=300)
plt.show()

# Dioscorea alata graph
plt.figure(figsize=(10, 5))
plt.bar(amino_acids, counts_protein_list, color="darkblue")

plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.title("Amino Acid Frequency: KAH7659748.1 t-snare proteins protein [Dioscorea alata]")

plt.tight_layout()
plt.savefig("dioscorea_amino_acid_frequency.png", dpi=300)
plt.show()

# print the counts
print("Translated protein counts:")
print(counts_translated)

print("\nProtein sequence counts:")
print(counts_protein)
