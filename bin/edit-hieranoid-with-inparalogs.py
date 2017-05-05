from Bio import SeqIO
from re import findall, match
from os import system, path, mkdir
from shutil import rmtree, copyfile

hieranoid = '/home/fatemeh/EnTrI/results/hieranoid/hieranoid-result.txt'
clusterspath = '/home/fatemeh/EnTrI/results/hieranoid/clusters-with-inparalogs.txt'
clusters = list()
with open(hieranoid) as from_file:
    for line in from_file:
        genes = findall('[\(,]([a-zA-Z0-9_]+):', line)
        clusters.append(genes)
with open(clusterspath, 'w') as clustersfile:
    for item in clusters:
        for ii in item:
            clustersfile.write(ii+'\t')
        clustersfile.write('\n')
