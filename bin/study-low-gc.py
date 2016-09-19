from Bio import SeqIO
from re import match, findall
from os import listdir, remove
from collections import defaultdict
from numpy import percentile


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../data/fasta-protein/chromosome/seqdb.fasta'
seqdict = read_fasta_sequences(seqdb)
biases = '../results/insertion-indices/check-biases'
lowgc = '../results/low-gc.txt'
highgc = '../results/high-gc.txt'
lowgcdesc = '../results/low-gc-with-description.txt'
highgcdesc = '../results/high-gc-with-description.txt'
gcinfo = '../results/low-gc-info.txt'
# highinfo = '../results/high-gc-info.txt'
list_of_files = listdir(biases)
with open(lowgc, 'w') as lowfile:
    with open(highgc, 'w') as highfile:
        for filename in list_of_files:
            gc = []
            with open(biases + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    gc.append(float(cells[4]))
            with open(biases+'/'+filename, 'r') as fromfile:
                for line in fromfile:
                    cells = line.split()
                    if float(cells[4]) <= percentile(gc, 25):
                        lowfile.write(cells[0]+'\t'+cells[5]+'\n')
                    elif percentile(gc, 25) < float(cells[4]) < percentile(gc, 75):
                        highfile.write(cells[0]+'\t'+cells[5]+'\n')
with open(lowgc, 'r') as fromfile:
    with open(lowgcdesc, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
            tofile.write(cells[0]+'\t'+cells[1]+'\t'+match_result.group(1)+'\n')
remove(lowgc)

with open(highgc, 'r') as fromfile:
    with open(highgcdesc, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]", seqdict[cells[0]].description)
            tofile.write(cells[0]+'\t'+cells[1]+'\t'+match_result.group(1)+'\n')
remove(highgc)

worddict = defaultdict(list)

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

with open(highgcdesc, 'r') as fromfile:
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

with open(gcinfo, 'w') as tofile:
    tofile.write('WORD\tLowGC\tHighGC\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))



# with open(lowinfo, 'w') as tofile:
#     with open(lowgcdesc, 'r') as fromfile:
#         for line in fromfile:
#             cells = line.split('\t')
#             description = cells[2].replace(',', ' ')
#             description = description.replace(':', ' ')
#             description = description.replace('_', ' ')
#             description = description.replace(';', ' ')
#             description = description.replace('-', ' ')
#             description = description.replace('(', ' ')
#             description = description.replace(')', ' ')
#             description = description.replace('.', ' ')
#             description = description.replace('/', ' ')
#             description = description.replace('%', ' ')
#             description = description.replace('*', ' ')
#             description = description.replace('\\', ' ')
#             description = description.replace('\'', ' ')
#             words = findall("\S+", description)
#             for item in words:
#                 if len(item) > 3:
#                     word = item.lower()
#                     if word in worddict.keys():
#                         worddict[word] += 1
#                     else:
#                         worddict[word] = 1
#     for item in sorted(worddict, key=worddict.get, reverse=True):
#         tofile.write(item+'\t'+str(worddict[item])+'\n')
#     tofile.write('SUM\t' + str(sum(worddict.values())))
#
# worddict = {}
# with open(highinfo, 'w') as tofile:
#     with open(highgcdesc, 'r') as fromfile:
#         for line in fromfile:
#             cells = line.split('\t')
#             description = cells[2].replace(',', ' ')
#             description = description.replace(':', ' ')
#             description = description.replace('_', ' ')
#             description = description.replace(';', ' ')
#             description = description.replace('-', ' ')
#             description = description.replace('(', ' ')
#             description = description.replace(')', ' ')
#             description = description.replace('.', ' ')
#             description = description.replace('/', ' ')
#             description = description.replace('%', ' ')
#             description = description.replace('*', ' ')
#             description = description.replace('\\', ' ')
#             description = description.replace('\'', ' ')
#             words = findall("\S+", description)
#             for item in words:
#                 if len(item) > 3:
#                     word = item.lower()
#                     if word in worddict.keys():
#                         worddict[word] += 1
#                     else:
#                         worddict[word] = 1
#     for item in sorted(worddict, key=worddict.get, reverse=True):
#         tofile.write(item+'\t'+str(worddict[item])+'\n')
#     tofile.write('SUM\t' + str(sum(worddict.values())))
