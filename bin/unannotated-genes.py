from re import match
from os import listdir, mkdir
from Bio import SeqIO


def read_fasta_sequences(seq):
    sequences = {}
    with open(seq, 'rU') as fasta_file:
            sequences.update(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')))
    return sequences

fastas = '../data/fasta-protein/chromosome/seqdb.fasta'
essentialities = '../results/biases/normalised-pca/'
newannot = '../results/unannotated/'
mkdir(newannot)
sequences = read_fasta_sequences(fastas)

list_of_files = listdir(essentialities)
for filename in list_of_files:
    with open(essentialities + filename, 'r') as fromfile:
        with open(newannot + filename, 'w') as tofile:
            tofile.write('ID\tname\tfunction\tcoord\tNPEQ\tessentiality\n')
            for line in fromfile:
                cells = line.split()
                if cells[0].endswith('added'):
                    desc = sequences[cells[0]].description
                    print(desc)
                    match_result = match(
                        '[^_]+_\d+added\s+\[[^/]+/(\S+)[^\]]+\]+\s+\[(.*)\]\s+\[.*\]\s+\[(.*)\]\s+\[.*\]', desc)
                    coord = match_result.group(1)
                    func = match_result.group(2)
                    name = match_result.group(3)
                    tofile.write(cells[0] + '\t' + name + '\t' + func + '\t' + coord + '\t' + cells[1] + '\t' +
                                 cells[2] + '\n')
