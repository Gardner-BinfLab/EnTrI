from Bio import SeqIO
from re import match, findall
from collections import defaultdict
from os import listdir


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

coreintersectpath = '../results/define-core-accessory-hieranoid/all/core-essential-genomes'
coreancestraliipath = '../results/define-core-accessory-hieranoid-ancestralii/all/core-essential-genomes'
coredollopath = '../results/define-core-accessory-hieranoid-fdollop/all/core-essential-genomes'
genomepath = '../data/fasta-protein/chomosome/seqdb.fasta'

ancintdes = '../results/word-enrichment/ancestral-intersection-description.txt'
nancintdes = '../results/word-enrichment/not-ancestral-intersection-description.txt'
ancintinfo = '../results/word-enrichment/ancestral-intersection-info.txt'

dolancdes = '../results/word-enrichment/dollo-ancestral-description.txt'
ndolancdes = '../results/word-enrichment/not-dollo-ancestral-description.txt'
dolancinfo = '../results/word-enrichment/dollo-ancestral-info.txt'

dolintdes = '../results/word-enrichment/dollo-intersection-description.txt'
ndolintdes = '../results/word-enrichment/not-dollo-intersection-description.txt'
dolintinfo = '../results/word-enrichment/dollo-intersection-info.txt'

list_of_files = listdir(coreintersectpath)
coreintersect = {}
for filename in list_of_files:
    if filename != "b.fasta":
        coreintersect.update(read_fasta_sequences(coreintersectpath+'/'+filename))

list_of_files = listdir(coreancestraliipath)
coreancestralii = {}
for filename in list_of_files:
    if filename != "b.fasta":
        coreancestralii.update(read_fasta_sequences(coreancestraliipath+'/'+filename))

list_of_files = listdir(coredollopath)
coredollo = {}
for filename in list_of_files:
    if filename != "b.fasta":
        coredollo.update(read_fasta_sequences(coredollopath+'/'+filename))

with open(nancintdes, 'w') as nfile:
    with open(ancintdes, 'w') as yfile:
        for item in coreancestralii.keys():
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", coreancestralii[item].description)
            if match_result:
                if item in coreintersect.keys():
                    yfile.write(item+'\t'+match_result.group(1)+'\n')
                else:
                    nfile.write(item+'\t'+match_result.group(1)+'\n')

worddict = defaultdict(list)

with open(ancintdes, 'r') as fromfile:
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

with open(nancintdes, 'r') as fromfile:
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

with open(ancintinfo, 'w') as tofile:
    tofile.write('WORD\tIntersection\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(ndolancdes, 'w') as nfile:
    with open(dolancdes, 'w') as yfile:
        for item in coredollo.keys():
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", coredollo[item].description)
            if match_result:
                if item in coreancestralii.keys():
                    yfile.write(item+'\t'+match_result.group(1)+'\n')
                else:
                    nfile.write(item+'\t'+match_result.group(1)+'\n')

worddict = defaultdict(list)

with open(dolancdes, 'r') as fromfile:
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

with open(ndolancdes, 'r') as fromfile:
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

with open(dolancinfo, 'w') as tofile:
    tofile.write('WORD\tAncestralii\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(ndolintdes, 'w') as nfile:
    with open(dolintdes, 'w') as yfile:
        for item in coredollo.keys():
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", coredollo[item].description)
            if match_result:
                if item in coreintersect.keys():
                    yfile.write(item+'\t'+match_result.group(1)+'\n')
                else:
                    nfile.write(item+'\t'+match_result.group(1)+'\n')

worddict = defaultdict(list)

with open(dolintdes, 'r') as fromfile:
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

with open(ndolintdes, 'r') as fromfile:
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

with open(dolintinfo, 'w') as tofile:
    tofile.write('WORD\tIntersection\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))