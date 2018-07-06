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

seqdb = '/home/fatemeh/EnTrI/results/deg/fasta-proteins/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/deg/clusters.txt'
outdir = '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid-fitch'
makedir(outdir)
speciestreedir = '/home/fatemeh/EnTrI/results/deg/trees/'
sequences = read_fasta_sequences(seqdb)
species_names = {'all': ['DEG1027', 'DEG1040', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1020', 'DEG1046', 'DEG1041',
                         'DEG1045', 'DEG1024', 'DEG1035', 'DEG1012', 'DEG1029', 'DEG1003', 'DEG1005', 'DEG1048',
                         'DEG1019', 'DEG1032', 'DEG1033', 'DEG1016', 'DEG1043', 'DEG1036', 'DEG1008', 'DEG1006',
                         'DEG1014', 'DEG1042', 'DEG1037', 'DEG1038', 'DEG1017', 'DEG1002', 'DEG1001'],
                 'proteobacteria': ['DEG1020', 'DEG1046', 'DEG1041', 'DEG1045', 'DEG1024', 'DEG1035', 'DEG1012',
                                    'DEG1029', 'DEG1003', 'DEG1005', 'DEG1048', 'DEG1019', 'DEG1032', 'DEG1033',
                                    'DEG1016', 'DEG1043', 'DEG1036', 'DEG1008'],
                 'gammaproteobacteria': ['DEG1012', 'DEG1029', 'DEG1003', 'DEG1005', 'DEG1048', 'DEG1019', 'DEG1032',
                                         'DEG1033', 'DEG1016', 'DEG1043', 'DEG1036']}
# species_names = ['A359', 'bbp', 'BOBLI757', 'BUMPF009', 'CWO', 'Sant', 'A35E', 'BCc', 'BPEN', 'BUMPG002', 'CWQ', 'SG',
#                  'AB162', 'BCHRO640', 'BTURN675', 'BUMPUSDA', 'CWS', 'SOPEG', 'BA000003', 'BCI', 'BUAMB', 'BUMPW106',
#                  'CWU', 'WIGMOR', 'BA000021', 'BCTU', 'BUAP5A', 'BUsg', 'IM45', 'BAKON', 'Bfl', 'BUAPTUC7', 'BVAF',
#                  'MEPCIT']
for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes'
    makedir(esscoredir)
    nesscoredir = speciesdir + '/non-essential-genomes'
    makedir(nesscoredir)

    counter = 0

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: {0} for i in range(num_species)}
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
                elif g.startswith('exDEG'):
                    s = g[2:9]
                else:
                    s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    if g.startswith('exDEG'):
                        gene_dict[s] = {0, 1}
                    else:
                        gene_dict[s] = {1}

            with open(speciestreedir + item + '.tre', 'r') as tree_file:
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

    gene_dict = {species_names[item][i]: {0} for i in range(num_species)}
    # ccc = 0
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
                elif g.startswith('exDEG'):
                    s = g[2:9]
                else:
                    s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    if g.startswith('exDEG'):
                        gene_dict[s] = {0, 1}
                    else:
                        gene_dict[s] = {1}
            if presence[index] == {1}:
                counter += 1
                with open(esscoredir + '/clust' + str(counter), 'w') as corefile:
                    for clustitem in list_of_genes:
                        corefile.write('>' + clustitem + '\n')
                        corefile.write(str(sequences[clustitem].seq) + '\n')

                # index += 1

            elif len(list_of_genes) > 0:
                counter += 1
                with open(nesscoredir + '/clust' + str(counter), 'w') as corefile:
                    for clustitem in list_of_genes:
                        corefile.write('>' + clustitem + '\n')
                        corefile.write(str(sequences[clustitem].seq) + '\n')
                        # if presence[index] == {0,1}:
                        #     ccc += 1
                        # index += 1
                        # print(ccc)
            index += 1