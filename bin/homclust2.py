# USAGE: homclust.py [-options] <seq_fasta_file> <seqdb_fasta_file> <output_directory>

# Needs Python, HMMER 3, and mafft installed.

# Clusters homologous proteins in different species. The first input file is a protein fasta file from one
# species and the second input file is a file resulted from merging all fasta files (including the first input file)
# for all given species.

# For merging fasta files, one can use:
# cat <file1> <file2> ... <filen> > seqdb.fasta

#! /bin/sh
""":"
exec python $0 ${1+"$@"}
"""

from os import path, mkdir, system, listdir
from shutil import rmtree, copytree
from Bio import SeqIO
from argparse import ArgumentParser
from math import sqrt
from timeit import default_timer
from dendropy import Tree


def parsetree(trepath):
    tree = Tree.get(path=trepath, schema='newick')
    distances = tree.phylogenetic_distance_matrix()  ##############https://pythonhosted.org/DendroPy/primer/phylogenetic_distances.html
    distance_matrix = list(list())
    for i, t1 in enumerate(tree.taxon_namespace):
        row = [None] * len(tree.taxon_namespace)
        for j, t2 in enumerate(tree.taxon_namespace):
            if i >= j:
                row[j] = float('inf')
            else:
                row[j] = distances(t1, t2)
        distance_matrix.append(row)
    nodes = tree.taxon_namespace.labels()
    return nodes, distance_matrix


def findclosest(matrix):
    minimum = float('inf')
    x, y = -1, -1
    for i, e1 in enumerate(matrix):
        for j, e2 in enumerate(e1):
            if matrix[i][j] < minimum:
                x = i
                y = j
                minimum = matrix[i][j]
    # print(indices)
    return x,y


def editdistances(names, matrix, sp):
    del matrix[sp]
    for row in matrix:
        del row[sp]
    del names[sp]
    return names, matrix


def mergefiles(file1, file2, merged):
    with open(merged, 'w') as tofile:
        write = tofile.write
        with open(file1, 'r') as fromfile:
            for line in fromfile:
                write(line)
    with open(file2, 'r') as fromfile:
        for line in fromfile:
            write(line)


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict


def extractfeatures(inpath, outpath):
    sequences = read_fasta_sequences(inpath)
    with open(outpath, 'w') as tofile:
        keys = sequences.keys()
        write = tofile.write
        for key in keys:
            write(sequences[key].id)
            seq = sequences[key].seq.upper()
            aadict = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0,
                      'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
            length = len(seq)
            for item in seq:
                if item in aadict.keys():
                    aadict[item] += 1
                elif item == 'B':
                    aadict['N'] += 0.5
                    aadict['D'] += 0.5
                elif item == 'Z':
                    aadict['Q'] += 0.5
                    aadict['E'] += 0.5
                elif item == 'J':
                    aadict['I'] += 0.5
                    aadict['L'] += 0.5
                else:
                    length -= 1
            sortedkeys = list(aadict.keys())
            sortedkeys.sort()
            for item in sortedkeys:
                if length == 0:
                    percentage = 0
                else:
                    percentage = float(aadict[item]) / length
                write('\t' + str(percentage))
            write('\n')


def readdata(inpath):
    names = list()
    data = list()
    with open(inpath, 'r') as fromfile:
        nappend = names.append
        dappend = data.append
        for line in fromfile:
            cells = line.split()
            nappend(cells[0])
            toappend = list(map(float, cells[1:]))
            dappend(toappend)
    return names, data


def initialclustering(inpath, distances, clusters):
    names, data = readdata(inpath)
    nrow = len(data)
    ncol = len(data[0])
    # distancematrix = sch.distance.pdist(data)
    # link = sch.linkage(distancematrix)
    # clusters = sch.fcluster(link,  200, criterion='maxclust')
    # return names, clusters
    with open(distances, 'w') as tofile:
        write = tofile.write
        for i in range(nrow):
            for j in range(i+1, nrow):
                # dist = 0
                # for k in range(ncol):
                #     dist += (data[i][k] - data[j][k])**2
                dist = sum([(data[i][k] - data[j][k])**2 for k in range(ncol)])
                dist = sqrt(dist)
                if dist < 0.05:
                    write(names[i] + '\t' + names[j] + '\t' + str(dist) + '\n')
    runmcl(distances, clusters)


def separatesequences(seqdb, clusters, outpath):
    sequences = read_fasta_sequences(seqdb)
    counter = 0
    with open(clusters, 'r') as fromfile:
        write = SeqIO.write
        for line in fromfile:
            cells = line.split()
            selectedsequences = [sequences[x] for x in cells]
            address = outpath + '/seq' + str(counter) + '.txt'
            write(selectedsequences, address, 'fasta')
            counter += 1


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


def runphmmer(eval, seqdb, outpath):
    try:
        system('phmmer -o /dev/null -E {0} --tblout {1} {2} {3}'.format(eval, outpath, seqdb, seqdb))
    except:
        raise SystemExit


def parsetblouts(inpath, outpath):
    edgelist = list()
    append = edgelist.append
    with open(outpath, 'w') as tofile:
        list_of_files = listdir(inpath)
        for filename in list_of_files:
            address = inpath + '/' + filename
            with open(address, 'r') as fromfile:
                for line in fromfile:
                    if not line.startswith('#'):
                        cells = line.split()
                        seq1 = min(cells[0], cells[2])
                        seq2 = max(cells[0], cells[2])
                        score = float(cells[5])
                        append([seq1, seq2, score])
        edgelist.sort()
        i = 0
        write = tofile.write
        while i < len(edgelist):
            seq1 = edgelist[i][0]
            seq2 = edgelist[i][1]
            if seq1 == seq2:
                i += 1
                continue
            score = float(edgelist[i][2])
            j = i + 1
            while j < len(edgelist) and edgelist[i][0:2] == edgelist[j][0:2]:
                score += float(edgelist[j][2])
                j += 1
            steps = j - i
            score /= steps
            i += steps
            write(seq1 + '\t' + seq2 + '\t' + str(score) + '\n')


def runmcl(inpath, outpath):
    try:
        system('mcl {0} --abc -o {1}'.format(inpath, outpath))
    except:
        raise SystemExit

start = default_timer()
et = 10**-10
parser = ArgumentParser(description='Clusters orthologous proteins. Needs python 3 or higher and HMMER 3.')
parser.add_argument('s', help='Directory containing all sequences in fasta format. The names of the files should be the'
                              'same as the name of the nodes in the species tree + \'.fa\'')
parser.add_argument('t', help='Rooted binary tree file in Newick format with non-zero branch lengths')
parser.add_argument('o', help='Path to the output')
parser.add_argument('-e', '--eva', help='E-value threshold for phmmer', default=et)
args = parser.parse_args()
seqpath = args.s
trepath = args.t
outpath = args.o
evalue = float(args.eva)

nodes, distancematrix = parsetree(trepath)
while len(distancematrix) > 1:
    sp1, sp2 = findclosest(distancematrix)
    seqdb = outpath + '/seq.txt'
    file1 = seqpath + nodes[sp1] + '.fa'
    file2 = seqpath + nodes[sp1] + '.fa'
    mergefiles(file1, file2, seqdb) #################correct the rest
    nodes, distancematrix= editdistances(nodes, distancematrix, sp2)
# seqdb = outpath + '/seq.txt'
# numspecies = mergefiles(seqpath, seqdb)
features = outpath + '/features.txt'
extractfeatures(seqdb, features)
distances = outpath + '/distances.txt'
initialclusters = outpath + '/initialclusters.txt'
initialclustering(features, distances, initialclusters)
sequences = outpath + '/sequences'
makedir(sequences)
separatesequences(seqdb, initialclusters, sequences)
phmmers = outpath + '/phmmers'
makedir(phmmers)
list_of_files = listdir(sequences)
for filename in list_of_files:
    runphmmer(eval, sequences + '/' + filename, phmmers + '/' + filename)
edges = outpath + '/edges.txt'
parsetblouts(phmmers, edges)
clusters = outpath + '/clusters.txt'
runmcl(edges, clusters)
stop = default_timer()
print(stop - start)