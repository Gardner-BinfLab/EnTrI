from Bio import SeqIO
from itertools import combinations
from re import search, findall, match
from os import listdir
from collections import defaultdict

def read_EC_files(filepath):
    ec_dict = defaultdict(list)
    locus_tags = {'t', 'ROD', 'ETEC', 'b', 'SL1344'}
    list_of_files = listdir(filepath)
    for filename in list_of_files:
        with open(filepath+filename) as ecfile:
            for line in ecfile:
                if line.startswith('ID'):
                    ec = []
                elif line.startswith('DE'):
                    ec += findall('\d+\.\d+\.\d+\.\d+', line)
                elif line.startswith('GN'):
                    search_results = search('OrderedLocusNames=([^;\{]*)', line)
                    if search_results:
                        find_results = findall('\w+', search_results.group(1))
                        for item1 in find_results:
                            for item2 in locus_tags:
                                if item1.startswith(item2):
                                    break
                                else:
                                    continue
                            break
                        if item1.startswith(item2):
                            ec_dict[item1] = ec
                            ec = []
    return ec_dict


def read_clusters(input_dir):
    sequence_occurrence = defaultdict(list)
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as cluster_file:
            for line in cluster_file:
                cells = line.split()
                sequence_occurrence[cells[1]].append(filename)
    return sequence_occurrence

ECpath = '../sequences/EC/'
clustdb = '../results/homclust/EFam-clusters'
ec_clusters = read_EC_files(ECpath)
efam_clusters = read_clusters(clustdb)
tp, fn, fp, tn = 0, 0, 0, 0
for key1, key2 in combinations(ec_clusters.keys(), r=2):
    if ec_clusters[key1] and ec_clusters[key2]:
        ecs1= []
        ecs2 = []
        for item in ec_clusters[key1]:
            match_results = match('\d+\.\d+\.\d+', item)
            ecs1.append(match_results.group(0))
        for item in ec_clusters[key2]:
            match_results = match('\d+\.\d+\.\d+', item)
            ecs2.append(match_results.group(0))
        ec_match = 0
        for item1 in ecs1:
            for item2 in ecs2:
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
            print '{}\t{}\t{}\t{}\t{}\t{}'.format(key1, ecs1, efam_clusters[key1], key2, ecs2, efam_clusters[key2])
        elif ec_match and not efam_match:
            fn += 1
        elif ec_match and efam_match:
            tp += 1
print tp, fn, fp, tn

# 8200 136616 608 6727354