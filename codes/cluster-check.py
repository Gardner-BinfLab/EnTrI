from Bio import SeqIO
from itertools import combinations
from re import search, findall
from os import listdir
from collections import defaultdict

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

def read_clusters(input_dir):
    sequence_occurrence = defaultdict(list)
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as cluster_file:
            for line in cluster_file:
                cells = line.split()
                sequence_occurrence[cells[1]].append(filename)
    return sequence_occurrence

seqdb = '../sequences/fasta-dna/chromosome/seqdb.fasta'
clustdb = '../results/homclust/EFam-clusters'
ec_clusters = read_fasta_sequences(seqdb)
efam_clusters = read_clusters(clustdb)
tp, fn, fp, tn = 0, 0, 0, 0
for key1, key2 in combinations(ec_clusters.keys(), r=2):
    search_results1 = search('\[EC_number=([^\]]*)\]', ec_clusters[key1].description)
    find_results1 = findall('(\d+\.\d+\.\d+)\.\d+', search_results1.group(1))
    search_results2 = search('\[EC_number=([^\]]*)\]', ec_clusters[key2].description)
    find_results2 = findall('(\d+\.\d+\.\d+)\.\d+', search_results2.group(1))
    if find_results1 and find_results2:
        ec_match = 0
        for item1 in find_results1:
            for item2 in find_results2:
                if item1 == item2:
                    ec_match = 1
                    break
            if ec_match:
                break
        efam_match = 0
        for item1 in efam_clusters[key1]:
            for item2 in efam_clusters[key2]:
                if item1 == item2:
                    efam_match = 1
                    break
            if efam_match:
                break
        if not ec_match and not efam_match:
            tn += 1
        elif not ec_match and efam_match:
            fp += 1
            print '{}\t{}\t{}\t{}\t{}\t{}'.format(key1, find_results1, efam_clusters[key1], key2, find_results2, efam_clusters[key2])
        elif ec_match and not efam_match:
            fn += 1
        elif ec_match and efam_match:
            tp += 1
print tp, fn, fp, tn

# 9138 158507 971 6838280
