#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match
from collections import defaultdict

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

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

seqdb = '/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/single-copy-genes'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-singcopy-ancestralii'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
list_of_files = listdir(clusters)
species_names = defaultdict()
species_names = {"all":["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "b"],"typhimurium":["STM", "SL1344", "STMMW"], "salmonella":["SEN", "SL1344", "STM", "STMMW", "t"],
    "ecoli":["CS17", "ETEC", "NCTC13441", "b"], "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "enterobacter":
    ["ENC"], "salmonellaecoli":["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441", "b"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441", "b", "ROD"],
    "klebsiellaenterobacter":["ERS227112", "BN373", "ENC"], "salmonellaty2":["t"], "salmonellap125109":["SEN"],
    "salmonellasl1344":["SL1344"], "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"],
    "ecolics17":["CS17"], "ecolih10407":["ETEC"], "ecolik12":["b"], "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"]}

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes'
    makedir(esscoredir)
    sesscoredir = speciesdir + '/core-sometimes-essential-genomes'
    makedir(sesscoredir)
    nesscoredir = speciesdir + '/core-never-essential-genomes'
    makedir(nesscoredir)
    essaccessorydir = speciesdir + '/accessory-essential-genomes'
    makedir(essaccessorydir)
    sessaccessorydir = speciesdir + '/accessory-sometimes-essential-genomes'
    makedir(sessaccessorydir)
    nessaccessorydir = speciesdir + '/accessory-never-essential-genomes'
    makedir(nessaccessorydir)

    esscoregenes = []
    nesscoregenes = []
    essaccessorygenes = []
    nessaccessorygenes = []

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: 0 for i in range(num_species)}
    essentiality_dict = {species_names[item][i]: 0 for i in range(num_species)}
    for filename in list_of_files:
        for key in gene_dict.keys():
            gene_dict[key] = 0
            essentiality_dict[key] = 0
        list_of_genes = []
        list_of_essential_genes = []

        with open('{0}/{1}'.format(clusters, filename)) as from_file:
            for line in from_file:
                cells = line.split()
                match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
                if not match_result:
                    match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
                if match_result:
                    name = match_result.group(1)
                    if name in gene_dict.keys():
                        list_of_genes.append(cells[1])
                        gene_dict[name] = 1
                        if 0 <= float(cells[4]):
                            essentiality_dict[name] = float(cells[4])
        list_of_species = []
        for species in gene_dict.keys():
            if gene_dict[species] == 1:
                list_of_species.append(species)
        print gene_dict
        ave = -1
        if len(essentiality_dict) > 0:
            if set(list_of_species) <= set(species_names["salmonella"]):
                ########################
            #####################
            #################
            #################
            ave = float(sum(essentiality_dict.values())) / len(essentiality_dict)
        if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
            if 0 <= ave <= 0.2:
                esscoregenes += list_of_genes
            else:
                nesscoregenes += list_of_genes
        else:
            if 0 <= ave <= 0.2:
                essaccessorygenes += list_of_genes
            else:
                nessaccessorygenes += list_of_genes
    esscoregenes = list(set(esscoregenes))
    esscoregenes.sort()
    nesscoregenes = list(set(nesscoregenes))
    nesscoregenes.sort()
    essaccessorygenes = list(set(essaccessorygenes))
    essaccessorygenes.sort()
    nessaccessorygenes = list(set(nessaccessorygenes))
    nessaccessorygenes.sort()

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


#fdollop -infile test.txt -intreefile test.treefile -outfile test.fdollop -method d -ancseq Y -treeprint N -progress N