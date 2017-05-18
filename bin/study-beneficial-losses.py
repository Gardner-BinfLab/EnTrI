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
essentiality = '../results/biases/dbscan'
loss = '../results/word-enrichment/beneficial-loss-description.txt'
nloss = '../results/word-enrichment/not-beneficial-loss-description.txt'
lossinfo = '../results/word-enrichment/beneficialloss-info.txt'

essen = '../results/word-enrichment/essential-description.txt'
nessen = '../results/word-enrichment/not-essential-description.txt'
esseninfo = '../results/word-enrichment/essential-info.txt'

none = '../results/word-enrichment/non-essential-description.txt'
nnone = '../results/word-enrichment/not-non-essential-description.txt'
noneinfo = '../results/word-enrichment/non-essential-info.txt'

list_of_files = listdir(essentiality)
with open(nloss, 'w') as nfile:
    with open(loss, 'w') as yfile:
        for filename in list_of_files:
            with open(essentiality + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
                    if match_result:
                        if cells[2] == 'beneficial-loss':
                            yfile.write(cells[0]+'\t'+match_result.group(1)+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + match_result.group(1) + '\n')

worddict = defaultdict(list)

with open(loss, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nloss, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
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




with open(nessen, 'w') as nfile:
    with open(essen, 'w') as yfile:
        for filename in list_of_files:
            with open(essentiality + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
                    if match_result:
                        if cells[2] == 'essential':
                            yfile.write(cells[0]+'\t'+match_result.group(1)+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + match_result.group(1) + '\n')

worddict = defaultdict(list)

with open(essen, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nessen, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(esseninfo, 'w') as tofile:
    tofile.write('WORD\tESSENTIAL\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))





with open(nnone, 'w') as nfile:
    with open(none, 'w') as yfile:
        for filename in list_of_files:
            with open(essentiality + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
                    if match_result:
                        if cells[2] == 'non-essential':
                            yfile.write(cells[0]+'\t'+match_result.group(1)+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + match_result.group(1) + '\n')

worddict = defaultdict(list)

with open(none, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nnone, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(noneinfo, 'w') as tofile:
    tofile.write('WORD\tNON_ESSENTIAL\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))