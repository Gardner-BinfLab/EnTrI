from itertools import combinations
from re import match, search
from os import listdir, mkdir

clusterspath = '../results/venn-entero-deg-endo/final-clusters.txt'
clusters_dict = {}
pathways = '../results/homologs.path'
eggnogs = '../results/homologs.eggnog'
groups = ['entero', 'endo', 'deg', 'prot', 'gamma']
for l in range(1,len(groups)+1):
    for subset in combinations(groups, l):
        subset = tuple(sorted(subset))
        clusters_dict[subset] = []

with open(clusterspath, 'r') as fromfile:
    for line in fromfile:
        line_dict = {item: 0 for item in groups}
        cells = line.split()
        for cell in cells[1:]:
            grp = match('([^-]*)', cell).group(1)
            line_dict[grp] += 1
        subset = [item for item in line_dict.keys() if line_dict[item] > 0]
        subset = tuple(sorted(subset))
        clusters_dict[subset].append(line.rstrip())

path_dict = {}
with open(pathways, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[2] = cells[2][1:-1]
        if cells[0].startswith('clust'):
            key = 'entero-'+cells[0]
            if key in path_dict.keys():
                path_dict[key] += ' /' + cells[2][0:-32]
            else:
                path_dict[key] = cells[2][0:-32]

eggnog_dict = {}
with open(eggnogs, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[1] = cells[1][1:-1]
        if len(cells[1]) == 0:
            cells[1] = '-'
        cells[2] = cells[2][1:-2]
        if len(cells[2]) == 0:
            cells[2] = '-'
        if cells[0].startswith('clust'):
            key = 'entero-' + cells[0]
            eggnog_dict[key] = cells[1] + '\t' + cells[2]

outpath = '../results/venn-entero-deg-endo/groups/'
mkdir(outpath)
for key in clusters_dict.keys():
    outpath = '../results/venn-entero-deg-endo/groups/'
    for item in key:
        outpath += item
    outpath += '.txt'
    if len(clusters_dict[key]) > 0:
        with open(outpath, 'w') as tofile:
            for item in clusters_dict[key]:
                if search('(entero-clust\d+)', item):
                    entero = search('(entero-clust\d+)', item).group(1)
                    if entero in path_dict.keys():
                        if entero in eggnog_dict.keys():
                            tofile.write(item + '\t' + path_dict[entero] + '\t' + eggnog_dict[entero] + '\n')
                        else:
                            tofile.write(item + '\t' + path_dict[entero] + '\t-\t-\n')
                    else:
                        if entero in eggnog_dict.keys():
                            tofile.write(item + '\t-' + '\t' + eggnog_dict[entero] + '\n')
                        else:
                            tofile.write(item + '\t-\t-\t-\n')
                else:
                    tofile.write(item + '\n')

