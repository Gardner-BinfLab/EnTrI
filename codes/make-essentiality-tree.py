#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path
from argparse import ArgumentParser
from re import match, compile

parser = ArgumentParser(description='Finds core genome and makes an essentiality tree from them')
parser.add_argument('indir', help='Directory of fasta clusters merged with essentiality information (merge-clust-plot/final_clusters)')
parser.add_argument('out', help='Output file name')
args = parser.parse_args()
input_dir = args.indir
out_path = args.out

# input_dir = '../results/merge-clust-plot/final_clusters'
# out_path = 'test.txt'
species_names = ["BN373", "CS17", "Ecoli9000q", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t"]
num_species = len(species_names)
gene_dict = {species_names[i]: 0 for i in range(num_species)}
ii_dict = {species_names[i]: i for i in range(num_species)}
name_dict = {"BN373": "KlePn.2", "CS17": "EscCo.1", "Ecoli9000q": "EscCo.2", "ENC": "EntCo.1",
             "ERS227112": "KlePn.1", "ETEC": "EscCo.3", "NCTC13441": "EscCo.4", "ROD": "CitRo.1",
             "SEN": "SalEn.1", "SL1344": "SalEn.3", "STM": "SalEn.2", "STMMW": "SalEn.4",
             "t": "SalEn.5", "b": "EscCo.5"}
iil = [None] * num_species

list_of_gene_files = []
list_of_insertion_indices = []

list_of_files = listdir(input_dir)
for filename in list_of_files:
    for key in gene_dict.keys():
        gene_dict[key] = 0
        iil[ii_dict[key]] = 0
    with open('{0}/{1}'.format(input_dir, filename)) as from_file:
        for line in from_file:
            cells = line.split()
            match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
            if not match_result:
                match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
            if match_result:
                name = match_result.group(1)
                if name in gene_dict.keys():
                    gene_dict[name] += 1
                    if float(cells[6]) < 0.2:
                        iil[ii_dict[name]] = 1
    #if gene_dict[min(gene_dict, key=gene_dict.get)] == 1 and gene_dict[max(gene_dict, key=gene_dict.get)] == 1 and 2 < sum(iil) < num_species - 1: # 2 <= sum(iil)-1 <= num_species - 3 # 1 < sum(iil)-1 < num_species-2 # 2 < sum(iil) < num_species-1
    if gene_dict[max(gene_dict, key=gene_dict.get)] == 1 and 2 < sum(iil) < num_species - 1:
        list_of_gene_files.append(path.basename(filename))
        list_of_insertion_indices.append(list(iil)) # Without list, it will append the reference to iil and if we change iil, the matrix will change.

num_genes = len(list_of_insertion_indices)
list_of_insertion_indices = [list(i) for i in zip(*list_of_insertion_indices)]
distance = [[0 for x in range(num_species)] for x in range(num_species)]
for i in range(num_species):
    for j in range(num_species):
        distance[i][j] = sum([list_of_insertion_indices[i][k]!=list_of_insertion_indices[j][k] for k in range(num_genes)]) # sum of xor

with open(out_path, 'w') as to_file:
    to_file.write('{0}\n'.format(num_species))
    for i in range(num_species):
        to_file.write('{0}   '.format(name_dict[species_names[i]]))
        for j in range(num_species - 1):
            to_file.write(str(distance[i][j])+'\t')
        to_file.write(str(distance[i][num_species - 1])+'\n')