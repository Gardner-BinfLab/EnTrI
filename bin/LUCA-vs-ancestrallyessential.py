from Bio import SeqIO
from os import listdir
from re import match, findall
from collections import defaultdict

#luca genes in gamma proteobacteria: 225

def read_annotations(filepath):
    annots = dict()
    with open(filepath, 'rU') as annotation_file:
        for line in annotation_file:
            cells = line.rstrip().split('\t')
            annots[cells[0]] = findall('\,(\w+)\@NOG', cells[8])[0] # in those genes with two cog IDs only the first one is in luca list
            # annots[cells[0]] = cells[-1]
    return annots

ancestralpath = '../results/luca-vs-ancestrallyessential/b.fasta.emapper.annotations'
lucapath = '../results/luca-vs-ancestrallyessential/luca-genes.txt'
allpath = '../results/eggnog-mapper/U00096.fasta.emapper.annotations'

ancestral = read_annotations(ancestralpath)
luca = []
all = read_annotations(allpath)
with open(lucapath, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('-'):
            luca.append(line.rstrip('\n'))

lucaancestral = {}
for key in ancestral.keys():
    if ancestral[key] in luca:
        lucaancestral[key] = ancestral[key]

print(lucaancestral)

# lucaall = {}
# for key in all.keys():
#     if all[key] in luca:
#         lucaall[key] = all[key]
#
# print(len(lucaall))