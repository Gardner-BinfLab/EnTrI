from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match, findall
from collections import defaultdict


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

def makedir(dirname):
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

seqdb = '../data/fasta-protein/chromosome/seqdb.fasta'
sequences = read_fasta_sequences(seqdb)
outdir = '../results/giant-tab/'
makedir(outdir)
clustdir = outdir + 'clusters/'
makedir(clustdir)
clusters = '../results/hieranoid/clusters.txt'
counter = 0
with open(clusters, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        cells = [i for i in cells if not i.startswith('ENC')]
        if len(cells) > 2:
            with open(clustdir + 'clust' + str(counter), 'w') as tofile:
                for gene in cells:
                    tofile.write('>' + gene + '\n')
                    tofile.write(str(sequences[gene].seq) + '\n')
                counter+=1