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
from shutil import rmtree
from Bio import SeqIO
from argparse import ArgumentParser
from random import randint
from re import findall
from timeit import default_timer


def addtostring(mainstring, position, addingstring):
    newstring = mainstring[:position] + addingstring + mainstring[position:]
    return newstring


def randomndigits(n):
    start = 10 ** (n - 1)
    end = (10 ** n) - 1
    return randint(start, end)


def addnodenames(tree):
    n = 10
    i = 0
    prevchar = ''
    while i < len(tree):
        char = tree[i]
        if char == ':' and prevchar == ')':
            nodename = 'n' + str(randomndigits(n))
            tree = addtostring(tree, i, nodename)
            i += len(nodename)
        elif char == ';' and prevchar == ')':
            nodename = 'n' + str(randomndigits(n)) + ':0'
            tree = addtostring(tree, i, nodename)
            i += len(nodename)
        prevchar = tree[i]
        i += 1
    return tree


def findnodes(tree):
    # findall_results = findall('[\(\),]([^\(^\)^,^:]*):([^,^\)^;]*)', tree)
    nodes = findall('[\(\),]([^\(^\)^,^:]*):', tree)
    return nodes


def calculaterootdistance(tree, nodes):
    rootdistance = {}
    for item in nodes:
        i = 0
        numpar = 0
        nodedistance = 0
        numparsaver = 0
        distance = ''
        nodename = ''
        char = tree[i]
        while i < len(tree):
            if char == '(':
                numpar += 1
                i += 1
                char = tree[i]
                while char != ':' and char != '(':
                    nodename += char
                    i += 1
                    char = tree[i]
            elif char == ')':
                if item == nodename or numpar < numparsaver:
                    nodedistance += float(distance)
                    numparsaver = numpar
                numpar -= 1
                i += 1
                char = tree[i]
                distance = ''
                nodename = ''
                while char != ':' and char != '(':
                    nodename += char
                    i += 1
                    char = tree[i]
            elif char == ':':
                i += 1
                char = tree[i]
                while char != ',' and char != ')' and char != ';':
                    distance += char
                    i += 1
                    char = tree[i]
            elif char == ',':
                if item == nodename or numpar < numparsaver:
                    nodedistance += float(distance)
                    numparsaver = numpar
                i += 1
                char = tree[i]
                distance = ''
                nodename = ''
                while char != ':' and char != '(':
                    nodename += char
                    i += 1
                    char = tree[i]
            elif char == ';':
                if item == nodename or numpar < numparsaver:
                    nodedistance += float(distance)
                    numparsaver = numpar
                i += 1
                distance = ''
                nodename = ''
            else:
                print('ERROR! The tree is not in a correct format.')
        rootdistance[item] = nodedistance
    return rootdistance


def calculatedistancematrix(tree, nodes, rootdistance):
    length = len(nodes)
    distancematrix = [[0 for x in range(length)] for y in range(length)]
    for x in range(length):
        for y in range(length):
            if tree.find(nodes[x]) <= tree.find(nodes[y]):
                item1 = nodes[x]
                item2 = nodes[y]
            else:
                item2 = nodes[x]
                item1 = nodes[y]
            i = 0
            numpar = 0
            numparsaver = 0
            nodename = ''
            char = tree[i]
            while i < len(tree):
                if char == '(':
                    numpar += 1
                    i += 1
                    char = tree[i]
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ')':
                    if item1 == nodename or numpar < numparsaver:
                        numparsaver = numpar
                    if item2 == nodename:
                        break
                    numpar -= 1
                    i += 1
                    char = tree[i]
                    nodename = ''
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ':':
                    i += 1
                    char = tree[i]
                    while char != ',' and char != ')' and char != ';':
                        i += 1
                        char = tree[i]
                elif char == ',':
                    if item1 == nodename or numpar < numparsaver:
                        numparsaver = numpar
                    if item2 == nodename:
                        break
                    i += 1
                    char = tree[i]
                    nodename = ''
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ';':
                    if item1 == nodename or numpar < numparsaver:
                        numparsaver = numpar
                    if item2 == nodename:
                        break
                    i += 1
                    nodename = ''

            while i < len(tree):
                if char == '(':
                    numpar += 1
                    i += 1
                    char = tree[i]
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ')':
                    if numpar < numparsaver:
                        break
                    numpar -= 1
                    i += 1
                    char = tree[i]
                    nodename = ''
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ':':
                    i += 1
                    char = tree[i]
                    while char != ',' and char != ')' and char != ';':
                        i += 1
                        char = tree[i]
                elif char == ',':
                    if numpar < numparsaver:
                        break
                    i += 1
                    char = tree[i]
                    nodename = ''
                    while char != ':' and char != '(':
                        nodename += char
                        i += 1
                        char = tree[i]
                elif char == ';':
                    if numpar < numparsaver:
                        break
                    i += 1
                    nodename = ''
            if item1 == item2:
                distancematrix[x][y] = 0.0
            elif nodename == '':
                distancematrix[x][y] = rootdistance[item1] - rootdistance[item2]
            else:
                distancematrix[x][y] = rootdistance[item1] + rootdistance[item2] - 2 * rootdistance[nodename]
    return distancematrix


def parsetree(inpath):
    with open(inpath, 'r') as fromfile:
        tree = fromfile.read()
    tree = addnodenames(tree)
    nodes = findnodes(tree)
    rootdistance = calculaterootdistance(tree, nodes)
    distancematrix = calculatedistancematrix(tree, nodes, rootdistance)
    print(nodes)
    for i in distancematrix:
        print(i)
    return nodes, distancematrix


def preprocess(matrix):
    maximum = float('inf')
    nrow = len(matrix)
    ncol = len(matrix[0])
    for i in range(nrow):
        for j in range(i, ncol):
            matrix[i][j] = maximum
    return matrix


def findindices(biglist, smalllist):
    indices = []
    for i in range(len(biglist)):
        if biglist[i] in smalllist:
            indices.append(i)
    return indices


def findmin(matrix, indices):
    minimum = float('inf')
    x, y = -1, -1
    for i in indices:
        for j in indices:
            if matrix[i][j] < minimum:
                x = i
                y = j
                minimum = matrix[i][j]
    return x,y


def findparent(biglist, smalllist, nodex, nodey, distancematrix):
    infinitive = float('inf')
    distancematrix[nodex][nodey] = infinitive



def findclosest(nodes, distancematrix, species):
    speciesind = findindices(nodes, species)
    minx, miny = findmin(distancematrix, speciesind)
    internals = list(set(nodes) - set(species))
    parent = findparent(nodes, internals, minx, miny, distancematrix)
    return nodes[minx], nodes[miny],



def mergefiles(inpath, outpath):
    list_of_files = listdir(inpath)
    with open(outpath, 'w') as tofile:
        write = tofile.write
        for filename in list_of_files:
            with open(inpath + '/' + filename, 'r') as fromfile:
                for line in fromfile:
                    write(line)
    numspecies = len(list_of_files)
    return numspecies


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
                dist **= 0.5
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


def runphmmer(evalue, seqdb, outpath):
    try:
        system('phmmer -o /dev/null -E {0} --tblout {1} {2} {3}'.format(evalue, outpath, seqdb, seqdb))
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
distancematrix = preprocess(distancematrix)
species = listdir(seqpath)
species = [item[:-3] for item in species]
sp1, sp2, parent = findclosest(nodes, distancematrix, species)
while minx != -1:
    distancematrix[minx][miny] = infinitive
    minx, miny = findmin(distancematrix, speciesind)
seqdb = outpath + '/seq.txt'
numspecies = mergefiles(seqpath, seqdb)
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
    runphmmer(evalue, sequences + '/' + filename, phmmers + '/' + filename)
edges = outpath + '/edges.txt'
parsetblouts(phmmers, edges)
clusters = outpath + '/clusters.txt'
runmcl(edges, clusters)
stop = default_timer()
print(stop - start)