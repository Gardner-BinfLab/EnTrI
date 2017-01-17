from os import path, mkdir, system, listdir
from shutil import rmtree
from Bio import SeqIO
from argparse import ArgumentParser
from timeit import default_timer
from dendropy import Tree
from itertools import combinations
from scipy.spatial.distance import pdist
import pandas as pd
import random
import string
import porc


def parse_tree(treepath, fmt):
    tree = Tree.get(path=treepath, schema=fmt)
    distances = tree.phylogenetic_distance_matrix()  ##############https://pythonhosted.org/DendroPy/primer/phylogenetic_distances.html
    nodes = tree.taxon_namespace.labels()
    distance_matrix = pd.DataFrame(float('inf'), index=nodes, columns=nodes)
    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for j, t2 in enumerate(tree.taxon_namespace[i+1:]):
            distance_matrix.at[t1.label, t2.label] = distances(t1, t2)
    return distance_matrix


def find_closest(distances):
    sp2 = distances.min().idxmin()
    sp1 = distances.idxmin()[sp2]
    return sp1, sp2


def edit_distances(distances, sp1, sp2):
    distances.ix[sp1 + ' ' + sp2, :] = distances.ix[[sp1, sp2], :].min(axis=0)
    distances.ix[:, sp1 + ' ' + sp2] = distances.ix[:, [sp1, sp2]].min(axis=1)
    distances.ix[sp1 + ' ' + sp2, sp1 + ' ' + sp2] = float('inf')
    distances.drop(sp1, axis=0, inplace=True)
    distances.drop(sp1, axis=1, inplace=True)
    distances.drop(sp2, axis=0, inplace=True)
    distances.drop(sp2, axis=1, inplace=True)
    return distances


# def add_sequence(sequences, outpath, counter, rand):
#     directory, filename = path.split(outpath)
#     mafftpath = directory + '/mafft' + rand + '_' + counter + '.msa'
#     print(mafftpath)
#     hmmbuildpath = directory + '/cluster' + rand + '_' + counter + '.hmm'
#     print(hmmbuildpath)
#     num_lines = sum(1 for line in open(sequences))
#     if num_lines > 2:
#         try:
#             system('mafft --text --quiet ' + sequences + ' > ' + mafftpath)
#             system('hmmbuild -o /dev/null --amino ' + hmmbuildpath + ' ' + mafftpath)
#             system('hmmemit -c ' + hmmbuildpath + ' >> ' + outpath)
#         except:
#             raise SystemExit
#     else:
#         with open(outpath, 'a') as tofile:
#             with open(sequences, 'r') as fromfile:
#                 for line in fromfile:
#                     tofile.write(line)
#
#
# def merge_sequences(clusters, sequences, outpath):
#     directory, filename = path.split(outpath)
#     genespath = directory + '/fastaseqs.fa'
#     rand = ''.join(random.choice(string.digits) for _ in range(6))
#     with open(clusters, 'r') as fromfile:
#         counter = 1
#         for line in fromfile:
#             with open(genespath, 'w') as tofile:
#                 write = tofile.write
#                 cells = line.split()
#                 for item in cells:
#                     write('>' + item + '\n' + str(sequences[item].seq) + '\n')
#             add_sequence(genespath, outpath, str(counter), rand)
#             counter += 1

# def sortout_clusters(temp, sp1, sp2, tempclusters, seq1, seq2):
#     sequence1 = read_fasta_sequences(seq1)
#     sequence2 = read_fasta_sequences(seq2)
#     clusters = temp + '/clusters_' + sp1 + '-' + sp2 + '.txt'
#     clust1 = temp + '/clusters_' + sp1 + '.txt'
#     clust2 = temp + '/clusters_' + sp2 + '.txt'
#     with open(tempclusters, 'r') as fromfile:
#         with open(clusters, 'w') as tofile:
#             for line in fromfile:
#                 cells = line.split()
#                 for item in cells:
#                     if item
#                 tofile.write('\n')


def orthology_call(seq, tree, outpath, thresh, evalue, schema, temp):
    distancematrix = parse_tree(tree, schema)
    temptemp = temp + '/tmp'
    tempfeatures = temp + '/features'
    templengths = temp + '/lengths'
    list_of_files = listdir(seq)
    sequences = {}
    for filename in list_of_files:
        sequences.update(porc.read_fasta_sequences(seq + '/' + filename))
    mkdir(tempfeatures)
    mkdir(templengths)
    while len(distancematrix) > 2:
        mkdir(temptemp)
        sp1, sp2 = find_closest(distancematrix)

        distancematrix = edit_distances(distancematrix, sp1, sp2)
        rmtree(temptemp)


def main():
    start = default_timer()
    parser = ArgumentParser(description='Clusters orthologous proteins. Needs python 3 or higher and HMMER 3.')
    parser.add_argument('s',
                        help='Directory containing all sequences in fasta format. The names of the files should be the'
                             'same as the name of the nodes in the species tree + \'.fa\'')
    parser.add_argument('p', help='Rooted binary phylogenetic tree file in Newick format with non-zero branch lengths')
    parser.add_argument('o', help='Path to the output')
    parser.add_argument('-f', '--fmt', help='Tree format, default = newick', default='newick')
    parser.add_argument('-e', '--eva', help='E-value threshold for phmmer', default=1e-10)
    parser.add_argument('-t', '--thr', help='Similarity threshold', default=0.05)
    args = parser.parse_args()
    seq = args.s
    tree = args.p
    outpath = args.o
    evalue = float(args.eva)
    thresh = float(args.thr)
    schema = args.fmt
    directory, filename = path.split(outpath)
    temp = directory + '/tmp' + ''.join(random.choice(string.digits) for _ in range(6))
    mkdir(temp)
    orthology_call(seq, tree, outpath, thresh, evalue, schema, temp)
    rmtree(temp)
    stop = default_timer()
    print(stop - start)


if __name__ == '__main__':
    main()