#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match, search
from collections import defaultdict


class Stack:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def size(self):
        return len(self.items)

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

def read_gene_essentiality(indir):
    list_of_files = listdir(indir)
    iidict = {}
    for filename in list_of_files:
        with open(indir+'/'+filename) as from_file:
            for line in from_file:
                cells = line.split()
                iidict[cells[0]] = cells[1]
    return iidict

def read_k12(inpath, iidict):
    with open(inpath) as from_file:
        for line in from_file:
            cells = line.split()
            iidict[cells[0]] = 'essential'
    return iidict

seqdb = '/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/hieranoid-result.txt'
insertion_indices = '/home/fatemeh/EnTrI/results/insertion-indices/gamma'
# k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid'
makedir(outdir)
speciestreedir = '/home/fatemeh/EnTrI/codes/speciestrees'
sequences = read_fasta_sequences(seqdb)
gene_essentiality = read_gene_essentiality(insertion_indices)
# gene_essentiality = read_k12(k12path, gene_essentiality)
species_names = defaultdict()
species_names = {"all":["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t"],"typhimurium":["STM", "SL1344", "STMMW"], "salmonella":["SEN", "SL1344", "STM", "STMMW", "t"],
    "ecoli":["CS17", "ETEC", "NCTC13441"], "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "enterobacter":
    ["ENC"], "salmonellaecoli":["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441", "ROD"],
    "klebsiellaenterobacter":["ERS227112", "BN373", "ENC"], "salmonellaty2":["t"], "salmonellap125109":["SEN"],
    "salmonellasl1344":["SL1344"], "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"],
    "ecolics17":["CS17"], "ecolih10407":["ETEC"], "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"]}

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes'
    makedir(esscoredir)
    nesscoredir = speciesdir + '/core-never-essential-genomes'
    makedir(nesscoredir)
    essaccessorydir = speciesdir + '/accessory-essential-genomes'
    makedir(essaccessorydir)
    nessaccessorydir = speciesdir + '/accessory-never-essential-genomes'
    makedir(nessaccessorydir)

    esscoregenes = []
    nesscoregenes = []
    essaccessorygenes = []
    nessaccessorygenes = []

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: 0 for i in range(num_species)}

    with open (clusters) as from_file:
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = 0
            list_of_genes = []
            stack_of_genes = Stack()
            sp_name = ''
            for char in line:
                if search('[(,]', char):
                    if sp_name in species_names[item]:

                    sp_name = ''
                elif search('[a-zA-Z0-9]', char):
                    sp_name += char
                elif search(')',char)



    esscoregenes = list(set(esscoregenes))
    esscoregenes.sort()
    nesscoregenes = list(set(nesscoregenes))
    nesscoregenes.sort()
    essaccessorygenes = list(set(essaccessorygenes))
    essaccessorygenes.sort()
    nessaccessorygenes = list(set(nessaccessorygenes))
    nessaccessorygenes.sort()

    prev_name = ''
    for gene in esscoregenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'esscoregenesfile' in locals():
                esscoregenesfile.close()
            esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
        esscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in nesscoregenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'nesscoregenesfile' in locals():
                nesscoregenesfile.close()
            nesscoregenesfile = open(nesscoredir + '/' + name + '.fasta', 'w')
        nesscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in essaccessorygenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'essaccessorygenesfile' in locals():
                essaccessorygenesfile.close()
            essaccessorygenesfile = open(essaccessorydir + '/' + name + '.fasta', 'w')
        essaccessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in nessaccessorygenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'nessaccessorygenesfile' in locals():
                nessaccessorygenesfile.close()
            nessaccessorygenesfile = open(nessaccessorydir + '/' + name + '.fasta', 'w')
        nessaccessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name