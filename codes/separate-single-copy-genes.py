#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from shutil import rmtree, copyfile
from re import match

def makedir(dirname):
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot'
outdir = '/home/fatemeh/EnTrI/results/single-copy-genes'
makedir(outdir)
species_names = ["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "b"]
num_species = len(species_names)
singcopy_clusters = []
list_of_files = listdir(clusters)
for filename in list_of_files:
    gene_dict = {species_names[i]: 0 for i in range(num_species)}
    with open('{0}/{1}'.format(clusters, filename)) as from_file:
        for line in from_file:
            cells = line.split()
            match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
            if not match_result:
                match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
            if match_result:
                name = match_result.group(1)
                if name in gene_dict.keys():
                    gene_dict[name] += 1
    if gene_dict[max(gene_dict, key=gene_dict.get)] == 1:
        singcopy_clusters.append(filename)
for filename in singcopy_clusters:
    copyfile(clusters + '/' + filename, outdir + '/' + filename)