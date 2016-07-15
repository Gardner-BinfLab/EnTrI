#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir, system
from Bio import SeqIO
from shutil import rmtree
from re import match, findall
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

def read_gene_essentiality(indir):
    list_of_files = listdir(indir)
    iidict = {}
    for filename in list_of_files:
        with open(indir+'/'+filename) as from_file:
            for line in from_file:
                cells = line.split()
                iidict[cells[0]] = cells[2]
    return iidict

def read_k12(inpath, iidict):
    with open(inpath) as from_file:
        for line in from_file:
            cells = line.split()
            iidict[cells[0]] = 'essential'
    return iidict

seqdb = '/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/hieranoid/hieranoid-result.txt'
insertion_indices = '/home/fatemeh/EnTrI/results/insertion-indices/normalised-insertion-indices'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fdollop'
speciestreedir = '/home/fatemeh/EnTrI/bin/speciestrees'
dollopout = '/home/fatemeh/EnTrI/bin/speciestrees/dollop.out'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
gene_essentiality = read_gene_essentiality(insertion_indices)
gene_essentiality = read_k12(k12path, gene_essentiality)
species_names = defaultdict()
species_names = {"all":["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "b", "SL3261"],"typhimurium":["STM", "SL1344", "STMMW", "SL3261"], "salmonella":["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"],
    "ecoli":["CS17", "ETEC", "NCTC13441", "b"], "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "enterobacter":
    ["ENC"], "salmonellacitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "ROD"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441", "b", "ROD"],
    "klebsiellaenterobacter":["ERS227112", "BN373", "ENC"], "salmonellaty2":["t"], "salmonellap125109":["SEN"],
    "salmonellasl1344":["SL1344"], "salmonellasl3261":["SL3261"], "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"],
    "ecolics17":["CS17"], "ecolih10407":["ETEC"], "ecolik12":["b"], "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"]}

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
    essentiality_dict = {species_names[item][i]: 0 for i in range(num_species)}
    nodeinfo = speciestreedir + '/nodeessentiality.txt'
    speciestree = speciestreedir + '/' + item + '.tre'

    with open (clusters) as from_file:
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = 0
                essentiality_dict[key] = 0
            list_of_genes = []

            # cells = line.split()
            findall_result = findall('(([a-zA-Z0-9]+?)_[a-zA-Z0-9]+):', line)
            temp_result = findall('(([a-zA-Z]+?)\d+):', line)
            for matches in temp_result:
                if matches[1] != 'n':
                    findall_result.append(matches)
            for element in findall_result:
                name = element[1]
                if name in gene_dict.keys():
                    list_of_genes.append(element[0])
                    gene_dict[name] = 1
                    if element[0] in gene_essentiality and gene_essentiality[element[0]] == 'essential':
                        essentiality_dict[name] = 1

            if len(species_names[item]) > 1:
                with open(nodeinfo, 'w') as nodefile:
                    nodefile.write('     ' + str(len(species_names[item])) + '    1\n')
                    for key in gene_dict:
                        if gene_dict[key] == 1:
                            nodefile.write(key + '          1\n')
                        else:
                            nodefile.write(key + '          0\n')

                try:
                    system(
                        'fdollop -infile {0} -intreefile {1} -outfile {2} -method d -progress N -treeprint N -ancseq Y'.format(
                            nodeinfo, speciestree, dollopout))
                except:
                    raise SystemExit
                presence = 0
                with open(dollopout, 'r') as dollopfile:
                    for line in dollopfile:
                        if line.startswith('root'):
                            cells = line.split()
                            if cells[2] == 'yes':
                                presence = 1

                with open(nodeinfo, 'w') as nodefile:
                    nodefile.write('     ' + str(len(species_names[item])) + '    1\n')
                    for key in essentiality_dict:
                        if essentiality_dict[key] == 1:
                            nodefile.write(key+'          1\n')
                        else:
                            nodefile.write(key+'          0\n')

                try:
                    system('fdollop -infile {0} -intreefile {1} -outfile {2} -method d -progress N -treeprint N -ancseq Y'.format(nodeinfo, speciestree, dollopout))
                except:
                    raise SystemExit
                essentiality = 0
                with open(dollopout, 'r') as dollopfile:
                    for line in dollopfile:
                        if line.startswith('root'):
                            cells = line.split()
                            if cells[2] == 'yes':
                                essentiality = 1

                if presence:
                    if essentiality:
                        esscoregenes += list_of_genes
                    else:
                        nesscoregenes += list_of_genes
                else:
                    if essentiality:
                        essaccessorygenes += list_of_genes
                    else:
                        nessaccessorygenes += list_of_genes
            else:
                name = list(gene_dict.keys())[0]
                if gene_dict[name]:
                    if essentiality_dict[name]:
                        esscoregenes += list_of_genes
                    else:
                        nesscoregenes += list_of_genes
                else:
                    if essentiality_dict[name]:
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
