from Bio import SeqIO


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '~/EnTrI/data/fasta-protein/chromosome/seqdb.fasta'
seqdict = read_fasta_sequences(seqdb)
lowgc = '~/EnTrI/results/low-gc.txt'
worddict = {}
with open(lowgc, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        