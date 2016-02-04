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
clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory'
# coredir = outdir + '/core-genomes'
esscoredir = outdir + '/core-essential-genomes'
sesscoredir = outdir + '/core-sometimes-essential-genomes'
nesscoredir = outdir + '/core-never-essential-genomes'
# accessorydir = outdir + '/accessory-genomes'
essaccessorydir = outdir + '/accessory-essential-genomes'
sessaccessorydir = outdir + '/accessory-sometimes-essential-genomes'
nessaccessorydir = outdir + '/accessory-never-essential-genomes'
makedir(outdir)
# makedir(coredir)
makedir(esscoredir)
makedir(sesscoredir)
makedir(nesscoredir)
# makedir(accessorydir)
makedir(essaccessorydir)
makedir(sessaccessorydir)
makedir(nessaccessorydir)
#core_genes = []
#core_essential_genes = []
esscoregenes = []
sesscoregenes = []
nesscoregenes = []
essaccessorygenes = []
sessaccessorygenes = []
nessaccessorygenes = []
species_names = ["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t", "b"]
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
                    if name in essentiality_dict.keys() and float(cells[4]) < 0.2:
                        essentiality_dict[name] = 1
    if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
        # core_genes += list_of_genes
        if essentiality_dict[min(essentiality_dict, key=essentiality_dict.get)] == 1:
            # core_essential_genes += list_of_genes
            esscoregenes += list_of_genes
        elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
            nesscoregenes += list_of_genes
        else:
            sesscoregenes += list_of_genes
    else:
        if gene_dict == essentiality_dict:
            essaccessorygenes += list_of_genes
        elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
            nessaccessorygenes += list_of_genes
        else:
            sessaccessorygenes += list_of_genes
# core_genes = list(set(core_genes))
# core_genes.sort()
# core_essential_genes = list(set(core_essential_genes))
# core_essential_genes.sort()
esscoregenes = list(set(esscoregenes))
esscoregenes.sort()
sesscoregenes = list(set(sesscoregenes))
sesscoregenes.sort()
nesscoregenes = list(set(nesscoregenes))
nesscoregenes.sort()
essaccessorygenes = list(set(essaccessorygenes))
essaccessorygenes.sort()
sessaccessorygenes = list(set(sessaccessorygenes))
sessaccessorygenes.sort()
nessaccessorygenes = list(set(nessaccessorygenes))
nessaccessorygenes.sort()
sequences = read_fasta_sequences(seqdb)

prev_name = ''
for item in esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'esscoregenesfile' in locals():
            esscoregenesfile.close()
        esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
    esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sesscoregenesfile' in locals():
            sesscoregenesfile.close()
        sesscoregenesfile = open(sesscoredir + '/' + name + '.fasta', 'w')
    sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'nesscoregenesfile' in locals():
            nesscoregenesfile.close()
        nesscoregenesfile = open(nesscoredir + '/' + name + '.fasta', 'w')
    nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'essaccessorygenesfile' in locals():
            essaccessorygenesfile.close()
        essaccessorygenesfile = open(essaccessorydir + '/' + name + '.fasta', 'w')
    essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sessaccessorygenesfile' in locals():
            sessaccessorygenesfile.close()
        sessaccessorygenesfile = open(sessaccessorydir + '/' + name + '.fasta', 'w')
    sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'nessaccessorygenesfile' in locals():
            nessaccessorygenesfile.close()
        nessaccessorygenesfile = open(nessaccessorydir + '/' + name + '.fasta', 'w')
    nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name