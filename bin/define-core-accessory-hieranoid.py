#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
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
clusters = '/home/fatemeh/EnTrI/results/hieranoid/clusters.txt'
insertion_indices = '/home/fatemeh/EnTrI/results/biases/dbscan'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
gene_essentiality = read_gene_essentiality(insertion_indices)
gene_essentiality = read_k12(k12path, gene_essentiality)
species_names = defaultdict()
species_names = {"all":["BN373", "ERS227112", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "b", "SL3261","BW25113", "EC958"],"typhimurium":["STM", "SL1344", "STMMW", "SL3261"],
    "salmonella":["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"], "ecoli":["NCTC13441", "b","BW25113", "EC958"],
    "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "salmonellacitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "ROD"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "NCTC13441", "b", "ROD","BW25113", "EC958"],
    "salmonellaty2":["t"], "salmonellap125109":["SEN"],
    "salmonellasl1344":["SL1344"], "salmonellasl3261":["SL3261"], "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"],
    "ecolik12":["b"], "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"],"ecoliBW25113":["BW25113"], "ecoliEC958":["EC958"]}

with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid/info.txt', 'w') as infofile:
    infofile.write('speciesname\tcoreessential\tcore\n')

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes/'
    makedir(esscoredir)
    nesscoredir = speciesdir + '/not-core-essential-genomes/'
    makedir(nesscoredir)

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: 0 for i in range(num_species)}
    essentiality_dict = {species_names[item][i]: 0 for i in range(num_species)}
    # for filename in list_of_files:
    with open (clusters) as from_file:
        counter = 0
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = 0
                essentiality_dict[key] = 0
            list_of_genes = []
            list_of_essential_genes = []

            genes = line.split()
            # findall_result = findall('(([a-zA-Z0-9]+?)_[a-zA-Z0-9]+):', line)
            # temp_result = findall('(([a-zA-Z]+?)\d+):', line)
            # for matches in temp_result:
            #     if matches[1] != 'n':
            #         findall_result.append(matches)
            for g in genes:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    gene_dict[s] = 1
                    if g in gene_essentiality and gene_essentiality[g] == 'essential':
                        list_of_essential_genes.append(g)
                        essentiality_dict[s] = 1
            # if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
            if sum(essentiality_dict.values()) >= len(gene_dict.keys()) * 80 / 100:
                with open(esscoredir + 'clust' + str(counter), 'w') as tofile:
                    for gene in list_of_genes:
                        tofile.write('>' + gene + '\n')
                        tofile.write(str(sequences[gene].seq) + '\n')
            else:
                with open(nesscoredir + 'clust' + str(counter), 'w') as tofile:
                    for gene in list_of_genes:
                        tofile.write('>' + gene + '\n')
                        tofile.write(str(sequences[gene].seq) + '\n')
            counter += 1