import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

df = pd.read_excel('02-Synechococcus elongatus PCC 7942.xlsx',sheet_name=0)
excluded = df.loc[df['essentiality'] == 'not_analyzed']['7942_ID']
# print(excluded.iloc[1])

fasta = read_fasta_sequences('../embl2peptides/02-Synechococcus-elongatus-PCC-7942.fasta')

with open('../excluded/02-Synechococcus-elongatus-PCC-7942.fasta', 'w') as tofile:
    for name in excluded:
        if name in fasta.keys():
            tofile.write(fasta[name].format('fasta'))