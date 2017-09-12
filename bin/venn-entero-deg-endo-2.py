from itertools import combinations
from re import match
from os import listdir

clusterspath = '../results/venn-entero-deg-endo/final-clusters.txt'
writepath = '../results/venn-entero-deg-endo/venn.txt'
singles_dict = {}
endopath = '../results/endosymbionts/define-core-accessory-hieranoid/core-essential-genomes/hhmake/'
files = listdir(endopath)
singles_dict['endo'] = len([i for i in files if i.startswith('endo-clust')])
enteropath = '../results/define-core-accessory-hieranoid-fitch-core-essentials/hhmake/'
files = listdir(enteropath)
singles_dict['entero'] = len([i for i in files if i.startswith('entero-clust')])
degpath = '../results/deg/define-core-accessory-hieranoid-fitch/all/core-essential-genomes/hhmake'
files = listdir(degpath)
singles_dict['deg'] = len([i for i in files if i.startswith('deg-clust')])
gammapath = '../results/deg/define-core-accessory-hieranoid-fitch/gammaproteobacteria/core-essential-genomes/hhmake'
files = listdir(gammapath)
singles_dict['gamma'] = len([i for i in files if i.startswith('gamma-clust')])
protpath = '../results/deg/define-core-accessory-hieranoid-fitch/proteobacteria/core-essential-genomes/hhmake'
files = listdir(protpath)
singles_dict['prot'] = len([i for i in files if i.startswith('prot-clust')])

clusters_dict = {}
groups = ['entero', 'endo', 'deg', 'prot', 'gamma']
for l in range(1,len(groups)+1):
    for subset in combinations(groups, l):
        subset = tuple(sorted(subset))
        clusters_dict[subset] = 0

line_dict = {item: 0 for item in groups}  # after each iteration over lines all values in this dict are zero
with open(clusterspath, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        for cell in cells[1:]:
            grp = match('([^-]*)', cell).group(1)
            line_dict[grp] += 1
            singles_dict[grp] -= 1
        morethan0 = [i for i in line_dict.keys() if line_dict[i] > 0]
        while len(morethan0) > 0:
            morethan0.sort()
            clusters_dict[tuple(morethan0)] += 1
            for item in morethan0:
                line_dict[item] -= 1
            morethan0 = [i for i in line_dict.keys() if line_dict[i] > 0]

for key in singles_dict.keys():
    clusters_dict[tuple([key])] += singles_dict[key]

with open(writepath, 'w') as tofile:
    for key in clusters_dict.keys():
        newkey = str(key)
        newkey = newkey[1:-1]
        newkey = newkey.replace('\'', '')
        newkey = newkey.replace(',', '')
        newkey = newkey.replace(' ', '&')
        newkey = '\'' + newkey + '\''
        newkey = newkey.replace('prot', 'Proteobacteria')
        newkey = newkey.replace('gamma', 'Gammaproteobacteria')
        newkey = newkey.replace('endo', 'Endosymbionts')
        newkey = newkey.replace('entero', 'Enterobacteriaceae')
        newkey = newkey.replace('deg', 'All')
        tofile.write(newkey + '\t' + str(clusters_dict[key]) + '\n')