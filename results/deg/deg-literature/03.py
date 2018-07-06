import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

df = pd.read_excel('03-Porphyromonas gingivalis strain ATCC 33277-short.xls',sheet_name=0)
excluded = df['Gene Locus']
# print(excluded.iloc[1])

fasta = read_fasta_sequences('../embl2peptides/03-Porphyromonas-gingivalis-strain-ATCC-33277.fasta')

with open('../excluded/03-Porphyromonas-gingivalis-strain-ATCC-33277.fasta', 'w') as tofile:
    for name in excluded:
        if name in fasta.keys():
            tofile.write(fasta[name].format('fasta'))