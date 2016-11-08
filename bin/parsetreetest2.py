from os import path, mkdir, system, listdir
from shutil import rmtree
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



trepath = '../../Desktop/testtt/tree.txt'
seqpath = '../../Desktop/testtt/genome2'
nodes, distancematrix = parsetree(trepath)
while len(distancematrix) > 1:
    sp1, sp2 = findclosest(distancematrix)
    nodes, distancematrix= editdistances(nodes, distancematrix, sp2)
# for row in distancematrix:
#     print(row)
# print(parent)