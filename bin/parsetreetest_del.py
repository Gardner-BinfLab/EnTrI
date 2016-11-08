from os import path, mkdir, system, listdir
from shutil import rmtree
#from Bio import SeqIO
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


# def calculatedistancematrix(tree, nodes, rootdistance):
#     length = len(nodes)
#     distancematrix = [[0 for x in range(length)] for y in range(length)]
#     for x in range(length):
#         for y in range(length):
#             item1 = nodes[min(x,y)]
#             item2 = nodes[max(x,y)]
#             i = 0
#             numpar = 0
#             numparsaver = 0
#             nodename = ''
#             char = tree[i]
#             while i < len(tree):
#                 if char == '(':
#                     numpar += 1
#                     i += 1
#                     char = tree[i]
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ')':
#                     if item1 == nodename or numpar < numparsaver:
#                         numparsaver = numpar
#                     if item2 == nodename:
#                         break
#                     numpar -= 1
#                     i += 1
#                     char = tree[i]
#                     nodename = ''
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ':':
#                     i += 1
#                     char = tree[i]
#                     while char != ',' and char != ')' and char != ';':
#                         i += 1
#                         char = tree[i]
#                 elif char == ',':
#                     if item1 == nodename or numpar < numparsaver:
#                         numparsaver = numpar
#                     if item2 == nodename:
#                         break
#                     i += 1
#                     char = tree[i]
#                     nodename = ''
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ';':
#                     if item1 == nodename or numpar < numparsaver:
#                         numparsaver = numpar
#                     if item2 == nodename:
#                         break
#                     i += 1
#                     nodename = ''
#
#             while i < len(tree):
#                 if char == '(':
#                     numpar += 1
#                     i += 1
#                     char = tree[i]
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ')':
#                     if numpar < numparsaver:
#                         break
#                     numpar -= 1
#                     i += 1
#                     char = tree[i]
#                     nodename = ''
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ':':
#                     i += 1
#                     char = tree[i]
#                     while char != ',' and char != ')' and char != ';':
#                         i += 1
#                         char = tree[i]
#                 elif char == ',':
#                     if numpar < numparsaver:
#                         break
#                     i += 1
#                     char = tree[i]
#                     nodename = ''
#                     while char != ':' and char != '(':
#                         nodename += char
#                         i += 1
#                         char = tree[i]
#                 elif char == ';':
#                     if numpar < numparsaver:
#                         break
#                     i += 1
#                     nodename = ''
#             if item1 == item2:
#                 distancematrix[x][y] = 0.0
#             elif nodename == '':
#                 distancematrix[x][y] = rootdistance[item1] - rootdistance[item2]
#             else:
#                 distancematrix[x][y] = rootdistance[item1] + rootdistance[item2] - 2 * rootdistance[nodename]
#     return distancematrix


def calculatedistancematrix(tree, nodes, rootdistance):
    length = len(nodes)
    distancematrix = [[float('inf') for x in range(length)] for y in range(length)]
    for x in range(length):
        for y in range(length):
            item1 = nodes[min(x,y)]
            item2 = nodes[max(x,y)]
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
    # print(indices)
    return x,y


def findparent(nodex, nodey, distancematrix):
    infinitive = float('inf')
    distancematrix[nodex][nodey] = infinitive
    distancematrix[nodey][nodex] = infinitive
    parentind1 = distancematrix[nodex].index(min(distancematrix[nodex]))
    column = [li[nodex] for li in distancematrix]
    parentind2 = column.index(min(column))
    if distancematrix[nodex][parentind1] < distancematrix[parentind2][nodex]:
        return parentind1
    else:
        return parentind2


def findclosest(nodes, distancematrix, species):
    speciesind = findindices(nodes, species)
    minx, miny = findmin(distancematrix, speciesind)
    return minx, miny


def editdistances(nodes, distancematrix, species, sp1, sp2, parent):
    # print(nodes)
    # print(nodes[sp1])
    # print(nodes[sp2])
    # for row in distancematrix:
    #     print(row)
    species.append(nodes[parent])
    species.remove(nodes[sp1])
    species.remove(nodes[sp2])
    del distancematrix[sp1]
    del distancematrix[sp2]
    for row in distancematrix:
        del row[sp1]
        del row[sp2]
    del nodes[sp1]
    del nodes[sp2]
    return nodes, distancematrix, species



trepath = '../../Desktop/testtt/tre2.txt'
seqpath = '../../Desktop/testtt/genome2'
nodes, distancematrix = parsetree(trepath)
distancematrix = preprocess(distancematrix)
for row in distancematrix:
    print(row)
print(nodes)
species = listdir(seqpath)
species = [item[:-3] for item in species]
while len(distancematrix) > 1:
    print(species)
    sp1, sp2 = findclosest(nodes, distancematrix, species)
    parent = findparent(sp1, sp2, distancematrix)
    nodes, distancematrix, species = editdistances(nodes, distancematrix, species, sp1, sp2, parent)
    # print(species)
# for row in distancematrix:
#     print(row)
# print(parent)