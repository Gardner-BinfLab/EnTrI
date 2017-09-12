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


seqdb = '/home/fatemeh/EnTrI/results/deg/fasta-proteins/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/deg/clusters.txt'
outdir = '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
species_names = {'all': ['DEG1001', 'DEG1008', 'DEG1020', 'DEG1029', 'DEG1038', 'DEG1045', 'DEG1002', 'DEG1012',
                         'DEG1021', 'DEG1031', 'DEG1039', 'DEG1046', 'DEG1003', 'DEG1013', 'DEG1023', 'DEG1034',
                         'DEG1040', 'DEG1047', 'DEG1005', 'DEG1014', 'DEG1024', 'DEG1035', 'DEG1041', 'DEG1006',
                         'DEG1015', 'DEG1027', 'DEG1036', 'DEG1042', 'DEG1007', 'DEG1017', 'DEG1028', 'DEG1037',
                         'DEG1043'], 'proteobacteria': ['DEG1035', 'DEG1024','DEG1012', 'DEG1005', 'DEG1003', 'DEG1029',
                                                        'DEG1015', 'DEG1036', 'DEG1013', 'DEG1043', 'DEG1020', 'DEG1046',
                                                        'DEG1041', 'DEG1045', 'DEG1028', 'DEG1031', 'DEG1008'],
                 'gammaproteobacteria': ['DEG1012', 'DEG1005', 'DEG1003', 'DEG1029', 'DEG1015', 'DEG1036', 'DEG1013',
                                         'DEG1043']}

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes/'
    makedir(esscoredir)
    sesscoredir = speciesdir + '/core-sometimes-essential-genomes/'
    makedir(sesscoredir)

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: 0 for i in range(num_species)}
    # for filename in list_of_files:
    with open (clusters) as from_file:
        counter = 1
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = 0
            list_of_genes = []

            genes = line.split()
            # findall_result = findall('(([a-zA-Z0-9]+?)_[a-zA-Z0-9]+):', line)
            # temp_result = findall('(([a-zA-Z]+?)\d+):', line)
            # for matches in temp_result:
            #     if matches[1] != 'n':
            #         findall_result.append(matches)
            for g in genes:
                s = g[0:7]
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    gene_dict[s] = 1
            # if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
            if sum(gene_dict.values()) >= len(gene_dict.keys())*80/100:
                with open(esscoredir + 'clust' + str(counter), 'w') as tofile:
                    for gene in list_of_genes:
                        tofile.write('>' + gene + '\n')
                        tofile.write(str(sequences[gene].seq) + '\n')
            else:
                with open(sesscoredir + 'clust' + str(counter), 'w') as tofile:
                    for gene in list_of_genes:
                        tofile.write('>' + gene + '\n')
                        tofile.write(str(sequences[gene].seq) + '\n')
            counter += 1
