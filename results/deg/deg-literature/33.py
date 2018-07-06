import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

excluded = ['Spy49_0071', 'Spy49_0210', 'Spy49_0235c', 'Spy49_0328', 'Spy49_0343', 'Spy49_0344', 'Spy49_0354c',
            'Spy49_0364c', 'Spy49_0407c', 'Spy49_0420', 'Spy49_0439', 'Spy49_0465', 'Spy49_0502c', 'Spy49_0565c',
            'Spy49_0786', 'Spy49_0798', 'Spy49_0847', 'Spy49_0917c', 'Spy49_1038c', 'Spy49_1042c', 'Spy49_1051c',
            'Spy49_1103', 'Spy49_1145c', 'Spy49_1217c', 'Spy49_1245', 'Spy49_1324c', 'Spy49_1463c', 'Spy49_1502c',
            'Spy49_1511c', 'Spy49_1516c', 'Spy49_1520c', 'Spy49_1528c', 'Spy49_1583c', 'Spy49_1643c', 'Spy49_1693',
            'Spy49_1705c', 'Spy49_1746c', 'Spy49_1771']

fasta = read_fasta_sequences('../embl2peptides/33-Streptococcus-pyogenes-NZ131.fasta')

with open('../excluded/33-Streptococcus-pyogenes-NZ131.fasta', 'w') as tofile:
    for name in excluded:
        if name in fasta.keys():
            tofile.write(fasta[name].format('fasta'))