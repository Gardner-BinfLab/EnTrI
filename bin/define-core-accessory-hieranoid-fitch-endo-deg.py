#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match


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
                if cells[2] == 'essential':
                        iidict[cells[0]] = {1}
                else:
                    iidict[cells[0]] = {0}
    return iidict

def read_k12(inpath, iidict):
    with open(inpath) as from_file:
        for line in from_file:
            cells = line.split()
            iidict[cells[0]] = {1}
    return iidict

seqdb = '/home/fatemeh/EnTrI/results/endosymbionts/fasta-proteins/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/endosymbionts/clusters.txt'
outdir = '/home/fatemeh/EnTrI/results/endosymbionts/define-core-accessory-hieranoid-fitch'
makedir(outdir)
speciestreedir = '/home/fatemeh/EnTrI/results/endosymbionts/tree.raxmlbootstrap'
sequences = read_fasta_sequences(seqdb)
counter = 0

esscoredir = outdir + '/core-essential-genomes'
makedir(esscoredir)
nesscoredir = outdir + '/non-essential-genomes'
makedir(nesscoredir)

esscoregenes = []
nesscoregenes = []
lenesscoregenes = 0
lennescoregenes = 0

# species_names = ['DEG1001', 'DEG1008', 'DEG1020', 'DEG1029', 'DEG1038', 'DEG1045', 'DEG1002', 'DEG1012', 'DEG1021',
#                  'DEG1031', 'DEG1039', 'DEG1046', 'DEG1003', 'DEG1013', 'DEG1023', 'DEG1034', 'DEG1040', 'DEG1047',
#                  'DEG1005', 'DEG1014', 'DEG1024', 'DEG1035', 'DEG1041', 'DEG1006', 'DEG1015', 'DEG1027', 'DEG1036',
#                  'DEG1042', 'DEG1007', 'DEG1017', 'DEG1028', 'DEG1037', 'DEG1043']
species_names = ['A359', 'bbp', 'BOBLI757', 'BUMPF009', 'CWO', 'Sant', 'A35E', 'BCc', 'BPEN', 'BUMPG002', 'CWQ', 'SG',
                 'AB162', 'BCHRO640', 'BTURN675', 'BUMPUSDA', 'CWS', 'SOPEG', 'BA000003', 'BCI', 'BUAMB', 'BUMPW106',
                 'CWU', 'WIGMOR', 'BA000021', 'BCTU', 'BUAP5A', 'BUsg', 'IM45', 'BAKON', 'Bfl', 'BUAPTUC7', 'BVAF',
                 'MEPCIT']
num_species = len(species_names)

gene_dict = {species_names[i]: {0} for i in range(num_species)}
presence = []

with open(clusters) as from_file:
    for line in from_file:
        for key in gene_dict.keys():
            gene_dict[key] = {0}
        list_of_genes = []
        stack_of_genes = Stack()

        genes = line.split()
        for g in genes:
            if g.startswith('DEG'):
                s = g[0:7]
            else:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
            if s in gene_dict.keys():
                list_of_genes.append(g)
                gene_dict[s] = {1}

        with open(speciestreedir , 'r') as tree_file:
            treeline = tree_file.readline()
            sp_name = ''
            i = 0
            char = treeline[i]
            while char != ';':
                if char == '(':
                    i += 1
                    char = treeline[i]
                elif char == ',':
                    if sp_name != '':
                        stack_of_genes.push(gene_dict[sp_name])
                        sp_name = ''
                    i += 1
                    char = treeline[i]
                elif char == ')':
                    if sp_name != '':
                        stack_of_genes.push(gene_dict[sp_name])
                        sp_name = ''

                    val1 = stack_of_genes.pop()
                    val2 = stack_of_genes.pop()
                    ave = val1 & val2
                    if len(ave) == 0:
                        ave = val1 | val2
                    stack_of_genes.push(ave)

                    i += 1
                    char = treeline[i]
                else:
                    sp_name += char
                    i += 1
                    char = treeline[i]

            presence.append(stack_of_genes.pop())

gene_dict = {species_names[i]: {0} for i in range(num_species)}
with open(clusters) as from_file:
    index = 0
    for line in from_file:
        for key in gene_dict.keys():
            gene_dict[key] = {0}
        list_of_genes = []
        genes = line.split()
        for g in genes:
            if g.startswith('DEG'):
                s = g[0:7]
            else:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
            if s in gene_dict.keys():
                list_of_genes.append(g)
                gene_dict[s] = {1}
        if presence[index] == {1}:
            counter += 1
            with open(esscoredir + '/clust' + str(counter), 'w') as corefile:
                for clustitem in list_of_genes:
                    corefile.write('>' + clustitem + '\n')
                    corefile.write(str(sequences[clustitem].seq) + '\n')

            index += 1

        else:
            counter += 1
            with open(nesscoredir + '/clust' + str(counter), 'w') as corefile:
                for clustitem in list_of_genes:
                    corefile.write('>' + clustitem + '\n')
                    corefile.write(str(sequences[clustitem].seq) + '\n')
            index += 1