import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

df = pd.read_excel('22-Acinetobacter baumannii.xlsx',sheet_name=0)
excluded = df.loc[pd.isnull(df['Input1 reads'])]['Gene accession no.']

fasta = read_fasta_sequences('../embl2peptides/22-Acinetobacter-baumannii.fasta')

with open('../excluded/22-Acinetobacter-baumannii.fasta', 'w') as tofile:
    for name in excluded:
        if name in fasta.keys():
            tofile.write(fasta[name].format('fasta'))