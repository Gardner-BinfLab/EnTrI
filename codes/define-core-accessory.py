#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

def makedir(dirname):
    # print ('Making directory \'{0}\'..'.format(dirname))
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

seqdb = '/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot/final_clusters'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory'
coredir = outdir + '/core-genomes'
esscoredir = outdir + '/core-essential-genomes'
accessorydir = outdir + '/accessory-genomes'
essaccessorydir = outdir + '/accessory-essential-genomes'
makedir(outdir)
makedir(coredir)
makedir(esscoredir)
mkdir(accessorydir)
mkdir(essaccessorydir)
core_genes = []
core_essential_genes = []
species_names = ["BN373", "CS17", "Ecoli9000q", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t", "b"]
num_species = len(species_names)
gene_dict = {species_names[i]: 0 for i in range(num_species)}
essentiality_dict = {species_names[i]: 0 for i in range(num_species-1)}

list_of_files = listdir(clusters)
for filename in list_of_files:
    for key in gene_dict.keys():
        gene_dict[key] = 0
        if key in essentiality_dict.keys():
            essentiality_dict[key] = 0
    list_of_genes = []
    with open('{0}/{1}'.format(clusters, filename)) as from_file:
        for line in from_file:
            cells = line.split()
            list_of_genes.append(cells[1])
            match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
            if not match_result:
                match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
            if match_result:
                name = match_result.group(1)
                if name in gene_dict.keys():
                    gene_dict[name] = 1
                    if name in essentiality_dict.keys() and float(cells[6]) < 0.2:
                        essentiality_dict[name] = 1
    if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
        core_genes += list_of_genes
        if essentiality_dict[min(essentiality_dict, key=essentiality_dict.get)] == 1:
            core_essential_genes += list_of_genes
core_genes = list(set(core_genes))
core_genes.sort()
core_essential_genes = list(set(core_essential_genes))
core_essential_genes.sort()
sequences = read_fasta_sequences(seqdb)

prev_name = ''
for item in core_genes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'coregenefile' in locals():
            coregenefile.close()
        coregenefile = open(coredir + '/' + name + '.fasta', 'w')
    coregenefile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in core_essential_genes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name in essentiality_dict.keys():
        if name != prev_name:
            if 'esscoregenefile' in locals():
                esscoregenefile.close()
            esscoregenefile = open(esscoredir + '/' + name + '.fasta', 'w')
        esscoregenefile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
all = list(sequences.keys())
all.sort()
for item in all:
    if item not in core_genes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', item)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', item)
        name = match_result.group(1)
        if name in gene_dict.keys():
            if name != prev_name:
                if 'accessoryfile' in locals():
                    accessoryfile.close()
                accessoryfile = open(accessorydir + '/' + name + '.fasta', 'w')
            accessoryfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in all:
    if item not in core_essential_genes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', item)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', item)
        name = match_result.group(1)
        if name in essentiality_dict.keys():
            if name != prev_name:
                if 'essaccessoryfile' in locals():
                    essaccessoryfile.close()
                essaccessoryfile = open(essaccessorydir + '/' + name + '.fasta', 'w')
            essaccessoryfile.write(sequences[item].format('fasta'))
    prev_name = name