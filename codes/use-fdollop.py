#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir, system
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
clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-dollo'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
list_of_files = listdir(clusters)
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
    sesscoregenes = []
    nesscoregenes = []
    essaccessorygenes = []
    sessaccessorygenes = []
    nessaccessorygenes = []

#fdollop -infile test.txt -intreefile test.treefile -outfile test.fdollop -method d -ancseq Y -treeprint N -progress N