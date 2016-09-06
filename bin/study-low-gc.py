from Bio import SeqIO
from re import match, findall
from os import listdir, remove


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../data/fasta-protein/chromosome/seqdb.fasta'
seqdict = read_fasta_sequences(seqdb)
biases = '../results/insertion-indices/check-biases'
lowgc = '../results/low-gc.txt'
lowgcdesc = '../results/low-gc-with-description.txt'
info = '../results/low-gc-info.txt'
list_of_files = listdir(biases)
with open(lowgc, 'w') as tofile:
    for filename in list_of_files:
        with open(biases+'/'+filename, 'r') as fromfile:
            for line in fromfile:
                cells = line.split()
                if float(cells[4]) <= 0.45:
                    tofile.write(cells[0]+'\t'+cells[5]+'\n')
with open(lowgc, 'r') as fromfile:
    with open(lowgcdesc, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
            tofile.write(cells[0]+'\t'+cells[1]+'\t'+match_result.group(1)+'\n')
remove(lowgc)

worddict = {}
with open(info, 'w') as tofile:
    with open(lowgcdesc, 'r') as fromfile:
        for line in fromfile:
            cells = line.split('\t')
            description = cells[2].replace(',', ' ')
            description = description.replace(':', ' ')
            description = description.replace('_', ' ')
            description = description.replace(';', ' ')
            description = description.replace('-', ' ')
            description = description.replace('(', ' ')
            description = description.replace(')', ' ')
            description = description.replace('.', ' ')
            description = description.replace('/', ' ')
            description = description.replace('%', ' ')
            description = description.replace('*', ' ')
            words = findall("\S+", description)
            for item in words:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word] += 1
                else:
                    worddict[word] = 1
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item+'\t'+str(worddict[item])+'\n')
    tofile.write('########## SUM: ' + str(sum(worddict.values())) + ' ##########')
