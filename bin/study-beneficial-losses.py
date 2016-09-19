from Bio import SeqIO
from os import listdir
from re import match, findall
from collections import defaultdict


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../data/fasta-protein/chromosome/seqdb.fasta'
seqdict = read_fasta_sequences(seqdb)
essentiality = '../results/insertion-indices/normalised-insertion-indices'
hloss = '../results/not-beneficial-loss-description.txt'
bloss = '../results/beneficial-loss-description.txt'
lossinfo = '../results/essentiality-info.txt'

list_of_files = listdir(essentiality)
with open(hloss, 'w') as hfile:
    with open(bloss, 'w') as bfile:
        for filename in list_of_files:
            with open(essentiality + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
                    if cells[2] == 'beneficial-loss':
                        bfile.write(cells[0]+'\t'+match_result.group(1)+'\n')
                    else:
                        hfile.write(cells[0] + '\t' + match_result.group(1) + '\n')

worddict = defaultdict(list)

with open(hloss, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
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
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(bloss, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
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
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(lossinfo, 'w') as tofile:
    tofile.write('WORD\tBeneficialLoss\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))
