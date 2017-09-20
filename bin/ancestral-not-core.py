from Bio import SeqIO
from re import match, findall
from collections import defaultdict
from os import listdir


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

ancestralpath = '../results/define-core-accessory-hieranoid-fitch/all/core-essential-genomes/b.fasta'
corepath = '../results/define-core-accessory-hieranoid-fasta/all/core-essential-genomes/b.fasta'
k12path = '../results/ecogene-k12.txt'
kegg = '../results/KEGG/escherichia_coli_K-12_MG1655.dat'
eggnog = '../results/eggnog-mapper/U00096.fasta.emapper.annotations'
writepath = '../results/core-ancestral/ancestral-not-core.tsv'
corewrite = '../results/core-ancestral/core.tsv'
ancestralwrite = '../results/core-ancestral/ancestral.tsv'
heatpath = '../results/interesting_genes/sometimes-essential-marked-dup-with-pathways.tsv'
heatwrite = '../results/core-ancestral/ancestral-not-core.heatmap'

ancestral = read_fasta_sequences(ancestralpath)
core = read_fasta_sequences(corepath)

k12 = []
with open(k12path, 'r') as fromfile:
    for line in fromfile:
        k12.append(line.strip())

path_dict = {}
with open(kegg, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[2] = cells[2][1:-1]
        if cells[0].startswith('b'):
            key = cells[0]
            if key in path_dict.keys():
                path_dict[key] += ' /' + cells[2][0:-32]
            else:
                path_dict[key] = cells[2][0:-32]

eggnog_dict = {}
with open(eggnog, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        if cells[0].startswith('b'):
            key = cells[0]
            eggnog_dict[key] = cells[4] + '\t' + cells[-1].strip()

with open (writepath, 'w') as tofile:
    tofile.write('gene_tag\tessentiality_in_K12\tKEGG_pathway\tgene_name\teggnog_annotation\n')
    for key in ancestral.keys():
        if key not in core.keys():
            if key in k12:
                ess = 'essential'
            else:
                ess = 'non-essential'
            if key in path_dict.keys() and key in eggnog_dict.keys():
                tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t' + eggnog_dict[key] + '\n')
            elif key in path_dict.keys():
                tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t-\t-' + '\n')
            elif key in eggnog_dict.keys():
                tofile.write(key + '\t' + ess + '\t-\t' + eggnog_dict[key] + '\n')
            else:
                tofile.write(key + '\t' + ess + '\t-\t-\t-' + '\n')

with open(ancestralwrite, 'w') as tofile:
    tofile.write('gene_tag\tessentiality_in_K12\tKEGG_pathway\tgene_name\teggnog_annotation\n')
    for key in ancestral.keys():
        if key in k12:
            ess = 'essential'
        else:
            ess = 'non-essential'
        if key in path_dict.keys() and key in eggnog_dict.keys():
            tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t' + eggnog_dict[key] + '\n')
        elif key in path_dict.keys():
            tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t-\t-' + '\n')
        elif key in eggnog_dict.keys():
            tofile.write(key + '\t' + ess + '\t-\t' + eggnog_dict[key] + '\n')
        else:
            tofile.write(key + '\t' + ess + '\t-\t-\t-' + '\n')

with open(corewrite, 'w') as tofile:
    tofile.write('gene_tag\tessentiality_in_K12\tKEGG_pathway\tgene_name\teggnog_annotation\n')
    for key in core.keys():
        if key in k12:
            ess = 'essential'
        else:
            ess = 'non-essential'
        if key in path_dict.keys() and key in eggnog_dict.keys():
            tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t' + eggnog_dict[key] + '\n')
        elif key in path_dict.keys():
            tofile.write(key + '\t' + ess + '\t' + path_dict[key] + '\t-\t-' + '\n')
        elif key in eggnog_dict.keys():
            tofile.write(key + '\t' + ess + '\t-\t' + eggnog_dict[key] + '\n')
        else:
            tofile.write(key + '\t' + ess + '\t-\t-\t-' + '\n')

with open(heatwrite, 'w') as tofile:
    with open(heatpath, 'r') as fromfile:
        line = fromfile.readline()
        tofile.write(line)
        for line in fromfile:
            cells = line.split('\t')
            if cells[29] in ancestral.keys() and cells[29] not in core.keys():
                tofile.write(line)