from os import path, mkdir, listdir
from subprocess import call
from shutil import rmtree, copyfile
from Bio import SeqIO
from argparse import ArgumentParser
from timeit import default_timer
from itertools import combinations, product, chain
from scipy.spatial.distance import pdist, cdist
import pandas as pd
# from numpy import average, std
# import matplotlib.pyplot as plt
import random
import string


def read_fasta_sequences(seqdir):
    list_of_files = listdir(seqdir)
    sequences = {}
    for filename in list_of_files:
        with open(seqdir + '/' + filename, 'rU') as fasta_file:
            sequences.update(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')))
    return sequences


def extract_features(sequences, chunksize, storepath):
    keys = list(sequences.keys())
    AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    cutoffs = list(range(0, len(keys), chunksize))
    maxlength = 30
    if len(keys) not in cutoffs:
        cutoffs.append(len(keys))
    with pd.get_store(storepath) as store:
        for i, s in enumerate(cutoffs[:-1]):
            e = cutoffs[i + 1]
            features = pd.DataFrame(0, index=keys[s:e], columns=AAs)
            lengths = pd.Series(0, index=keys[s:e])
            for key in keys[s:e]:
                seq = sequences[key].seq.upper()
                length = len(seq)
                for item in seq:
                    if item in AAs:
                        features.at[key, item] += 1
                    elif item == 'B':
                        features.at[key, 'N'] += 0.5
                        features.at[key, 'D'] += 0.5
                    elif item == 'Z':
                        features.at[key, 'Q'] += 0.5
                        features.at[key, 'E'] += 0.5
                    elif item == 'J':
                        features.at[key, 'I'] += 0.5
                        features.at[key, 'L'] += 0.5
                    else:
                        length -= 1
                lengths.at[key] = length
            features = features.divide(lengths, axis='index')
            store.append('features', features, min_itemsize={'index': maxlength})
            store.append('lengths', lengths, min_itemsize={'index': maxlength})


def calculate_distances(storepath, lenfeatures, chunksize, thresh, outpath):
    cutoffs = list(range(0, lenfeatures, chunksize))
    if lenfeatures not in cutoffs:
        cutoffs.append(lenfeatures)
    with pd.get_store(storepath) as store:
        for i, s1 in enumerate(cutoffs[:-1]):
            e1 = cutoffs[i + 1]
            for j, s2 in enumerate(cutoffs[i:-1]):
                e2 = cutoffs[i + j + 1]
                if s1 == s2:
                    chunk = store.select('features', start=s1, stop=e1)
                    distances = pdist(chunk, 'euclidean')
                    tuples = combinations(chunk.index, 2)
                else:
                    chunk1 = store.select('features', start=s1, stop=e1)
                    chunk2 = store.select('features', start=s2, stop=e2)
                    distances = cdist(chunk1, chunk2, 'euclidean')
                    distances = list(chain.from_iterable(distances))
                    tuples = product(chunk1.index, chunk2.index)
                tuplesdf = pd.DataFrame.from_records(tuples)
                distancesdf = pd.Series(distances)
                tupledist = pd.concat([tuplesdf, distancesdf], axis=1)
                tupledist.columns=['g1', 'g2', 'dist']
                tuplesmalldist = tupledist.loc[tupledist['dist'] < thresh]
                tuplesmalldist.to_csv(outpath, sep='\t', header=False, index=False, mode='a')


def run_mcl(inpath, outpath):
    try:
        call(['mcl', inpath, '--abc', '-o', outpath], shell=False)
    except:
        raise SystemExit('Unable to run MCL. Check if:\n 1- ' + inpath + ' exists.\n 2- MCL is installed.')


def find_centres(clusters, sequences, outpath):
    AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    with open(clusters, 'r') as fromfile:
        with pd.get_store(outpath) as store:
            for line in fromfile:
                keys = line.split()
                length = 0
                centres = pd.Series(0, index=AAs)
                for key in keys:
                    seq = sequences[key].seq.upper()
                    length += len(seq)
                    for item in seq:
                        if item in AAs:
                            centres.at[item] += 1
                        elif item == 'B':
                            centres.at['N'] += 0.5
                            centres.at['D'] += 0.5
                        elif item == 'Z':
                            centres.at['Q'] += 0.5
                            centres.at['E'] += 0.5
                        elif item == 'J':
                            centres.at['I'] += 0.5
                            centres.at['L'] += 0.5
                        else:
                            length -= 1
                centres = centres.divide(length, axis='index')
                store.append('centres', centres, min_itemsize={'index': maxlength})



def pairwise_orthology_call(seq, outpath, chunksize, thresh, temp):
    sequences = read_fasta_sequences(seq)
    storepath = temp + '/porc.h5'
    extract_features(sequences, chunksize, storepath)
    distances = temp + '/distances.txt'
    lenfeatures = len(list(sequences.keys()))
    calculate_distances(storepath, lenfeatures, chunksize, thresh, distances)
    initclusters = temp + '/' + 'initial_clusters.txt'
    run_mcl(distances, initclusters)
    find_centres(initclusters, sequences, outpath)



def main():
    start = default_timer()
    parser = ArgumentParser(description='Clusters orthologous proteins. Needs Python 2.7 or higher and MCL.')
    parser.add_argument('p', help='Directory containing all proteomes in fasta format')
    parser.add_argument('o', help='Path to the output')
    parser.add_argument('-s', '--thr', help='Similarity threshold', default=0.05)
    parser.add_argument('-c', '--chnk', help='Chunk size', default=5000)
    args = parser.parse_args()
    seq = args.p
    outpath = args.o
    thresh = float(args.thr)
    chunksize = int(args.chnk)
    directory, filename = path.split(outpath)
    temp = directory + '/tmp' + ''.join(random.choice(string.digits) for _ in range(6))
    mkdir(temp)
    pairwise_orthology_call(seq, outpath, chunksize, thresh, temp)
    #rmtree(temp)
    stop = default_timer()
    print(stop - start)


if __name__ == '__main__':
    main()###http://pandas-docs.github.io/pandas-docs-travis/io.html#io-hdf5
