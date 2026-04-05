# Biopython-Analysis
Sequence analysis pipeline using Biopython

This project analyzes and visualizes the amino acid composition of t-SNARE proteins from two plant species:

Landoltia punctata (translated from nucleotide sequence)
Dioscorea alata (protein sequence)

Using Biopython and Matplotlib, the script:
- Retrieves sequences from NCBI using Entrez
- Translates nucleotide sequences into proteins
- Counts amino acid frequencies
- Generates bar charts for visualization

Nucleotide sequence: JZ990160.1 (Landoltia punctata)
Protein sequence: KAH7659748.1 (Dioscorea alata)

Requirements: 
Install the required Python libraries
pip install biopython matplotlib

How To Run:
Clone this repository:
git clone https://github.com/sahanaarvind2010/Biopython-Analysis.git
Navigate into the folder:
cd your-repo-name
Run the script:
python amino_acid_frequency_t_snare.py

Output: 
The script generates two bar charts:

Landoltia punctata (translated protein)
File: landoltia_amino_acid_frequency.png
Dioscorea alata protein
File: dioscorea_amino_acid_frequency.png

Each chart shows the frequency of the 20 standard amino acids.

Author: 
Sahana Arvind 

This project is for educational and research purposes.
