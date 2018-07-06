import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

df = pd.read_excel('01-Mycobacterium tuberculosis H37Rv.XLSX',sheet_name=0)
excluded = df.loc[df['Call'] == 'short']['Rv Number']
# print(excluded.iloc[1])

fasta = read_fasta_sequences('../embl2peptides/01-Mycobacterium-tuberculosis-H37Rv.fasta')

with open('../excluded/01-Mycobacterium-tuberculosis-H37Rv.fasta', 'w') as tofile:
    for name in excluded:
        tofile.write(fasta[name].format('fasta'))