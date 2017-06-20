from Bio import SeqIO
from os import listdir
from re import match, findall
from collections import defaultdict


def read_annotations(filepath):
    annots = dict()
    with open(filepath, 'rU') as annotation_file:
        for line in annotation_file:
            cells = line.rstrip().split('\t')
            # annots[cells[0]] = 'ENOG41' + findall('\,(\w+)\@gproNOG', cells[8])[0]
            annots[cells[0]] = cells[-1]
    return annots

seqdb = '../results/eggnog-mapper/seqdb.fasta.emapper.annotations'
seqdict = read_annotations(seqdb)
essentiality = '../results/biases/dbscan'
loss = '../results/word-enrichment2/beneficial-loss-description.txt'
nloss = '../results/word-enrichment2/not-beneficial-loss-description.txt'
lossinfo = '../results/word-enrichment2/beneficialloss-info.txt'

essen = '../results/word-enrichment2/essential-description.txt'
nessen = '../results/word-enrichment2/not-essential-description.txt'
esseninfo = '../results/word-enrichment2/essential-info.txt'

none = '../results/word-enrichment2/non-essential-description.txt'
nnone = '../results/word-enrichment2/not-non-essential-description.txt'
noneinfo = '../results/word-enrichment2/non-essential-info.txt'

list_of_files = listdir(essentiality)
with open(nloss, 'w') as nfile:
    with open(loss, 'w') as yfile:
        for filename in list_of_files:
            with open(essentiality + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    # match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]])
                    # if match_result:
                    if cells[0] in seqdict.keys():
                        if cells[2] == 'beneficial-loss':
                            yfile.write(cells[0]+'\t'+seqdict[cells[0]]+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + seqdict[cells[0]] + '\n')

worddict = defaultdict(list)

with open(loss, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
        if word in worddict.keys():
            worddict[word][0] += 1
        else:
            worddict[word] = [1, 0]

with open(nloss, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
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
                    # match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]])
                    # if match_result:
                    if cells[0] in seqdict.keys():
                        if cells[2] == 'essential':
                            yfile.write(cells[0]+'\t'+seqdict[cells[0]]+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + seqdict[cells[0]] + '\n')

worddict = defaultdict(list)

with open(essen, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
        if word in worddict.keys():
            worddict[word][0] += 1
        else:
            worddict[word] = [1, 0]

with open(nessen, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
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
                    # match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]])
                    # if match_result:
                    if cells[0] in seqdict.keys():
                        if cells[2] == 'non-essential':
                            yfile.write(cells[0]+'\t'+seqdict[cells[0]]+'\n')
                        else:
                            nfile.write(cells[0] + '\t' + seqdict[cells[0]] + '\n')

worddict = defaultdict(list)

with open(none, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
        if word in worddict.keys():
            worddict[word][0] += 1
        else:
            worddict[word] = [1, 0]

with open(nnone, 'r') as fromfile:
    for line in fromfile:
        cells = line.rstrip().split('\t')
        word = cells[1]
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