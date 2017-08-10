from argparse import ArgumentParser
parser = ArgumentParser(description='Clusters orthologous proteins. Needs Python 2.7 or higher, MCL, and HMMER 3')
parser.add_argument('-i', help='File containing initial clusters')
parser.add_argument('-m', help='Calculate means of clusters', action='store_true')
args = parser.parse_args()
clust = args.i
means = args.m
print(means)