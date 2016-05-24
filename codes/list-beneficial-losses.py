from Bio import SeqIO
from os import listdir

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../sequences/fasta-protein/chromosome/seqdb.fasta'
iis = '../results/check-biases/without-ends.txt'
coreesspath = '../results/define-core-accessory/all/core-essential-genomes'
corenespath = '../results/define-core-accessory/all/core-never-essential-genomes'
accessoryesspath = '../results/define-core-accessory/all/accessory-essential-genomes'
accessorynespath = '../results/define-core-accessory/all/accessory-never-essential-genomes'
benloss = '../results/beneficial-losses.fasta'
coreess = '../results/core-essential.fasta'
corenes = '../results/core-never-essential.fasta'
accessoryess = '../results/accessory-essential.fasta'
accessorynes = '../results/accessory-never-essential.fasta'
sequence_dict = read_fasta_sequences(seqdb)
with open(iis, 'r') as iisfile:
    with open(benloss, 'w') as benlossfile:
        for line in iisfile:
            cells = line.split()
            if float(cells[1]) > 2:
                benlossfile.write('>')
                benlossfile.write(sequence_dict[cells[0]].description)
                benlossfile.write('\n')
                benlossfile.write(str(sequence_dict[cells[0]].seq))
                benlossfile.write('\n')


list_of_files = listdir(coreesspath)
with open(coreess, 'w') as outfile:
    for filename in list_of_files:
        with open(coreesspath+'/'+filename) as infile:
            for line in infile:
                outfile.write(line)

list_of_files = listdir(corenespath)
with open(corenes, 'w') as outfile:
    for filename in list_of_files:
        with open(corenespath+'/'+filename) as infile:
            for line in infile:
                outfile.write(line)

list_of_files = listdir(accessoryesspath)
with open(accessoryess, 'w') as outfile:
    for filename in list_of_files:
        with open(accessoryesspath+'/'+filename) as infile:
            for line in infile:
                outfile.write(line)

list_of_files = listdir(accessorynespath)
with open(accessorynes, 'w') as outfile:
    for filename in list_of_files:
        with open(accessorynespath+'/'+filename) as infile:
            for line in infile:
                outfile.write(line)