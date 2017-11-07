from itertools import combinations
from re import match
from os import listdir

writepath_init = '../results/venn-entero-deg-endo/final-clusters-plus-single-groups.txt'
writepath_final = '../results/venn-entero-deg-endo/final-clusters-plus-single-groups.fasta'
clusterspath = '../results/venn-entero-deg-endo/final-clusters.txt'
singles_dict = {}
endopath = '../results/endosymbionts/define-core-accessory-hieranoid/core-essential-genomes/hhmake/'
files = listdir(endopath)
singles_dict['endo'] = [i.split('.')[0] for i in files if i.startswith('endo-clust')]
enteropath = '../results/define-core-accessory-hieranoid/all/core-essential-genomes/hhmake/'
files = listdir(enteropath)
singles_dict['entero'] = [i.split('.')[0] for i in files if i.startswith('entero-clust')]
degpath = '../results/deg/define-core-accessory-hieranoid/all/core-essential-genomes/hhmake'
files = listdir(degpath)
singles_dict['deg'] = [i.split('.')[0] for i in files if i.startswith('deg-clust')]
gammapath = '../results/deg/define-core-accessory-hieranoid/gammaproteobacteria/core-essential-genomes/hhmake'
files = listdir(gammapath)
singles_dict['gamma'] = [i.split('.')[0] for i in files if i.startswith('gamma-clust')]
protpath = '../results/deg/define-core-accessory-hieranoid/proteobacteria/core-essential-genomes/hhmake'
files = listdir(protpath)
singles_dict['prot'] = [i.split('.')[0] for i in files if i.startswith('prot-clust')]

with open(clusterspath, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        for cell in cells[1:]:
            grp = match('([^-]*)', cell).group(1)
            singles_dict[grp].remove(cell)

print(singles_dict)

counter = 1
with open (writepath_init, 'w') as tofile:
    with open (clusterspath, 'r') as fromfile:
        for line in fromfile:
           tofile.write(line)
           counter += 1
    for key in singles_dict.keys():
        for item in singles_dict[key]:
            tofile.write('clust-' + str(counter) + '\t' + item + '\n')
            counter += 1

addresses = {'endo': '../results/endosymbionts/define-core-accessory-hieranoid/core-essential-genomes/',
             'prot': '../results/deg/define-core-accessory-hieranoid/proteobacteria/core-essential-genomes/',
             'gamma': '../results/deg/define-core-accessory-hieranoid/gammaproteobacteria/core-essential-genomes/',
             'entero': '../results/define-core-accessory-hieranoid/all/core-essential-genomes/',
             'deg': '../results/deg/define-core-accessory-hieranoid/all/core-essential-genomes/'}
with open(writepath_init, 'r') as fromfile:
    with open(writepath_final, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            tofile.write('>'+cells[0]+'\n')
            clust = cells[1].split('-')
            print(cells[1])
            with open(addresses[clust[0]] + clust[1], 'r') as clustfile:
                clustfile.readline()
                seq = clustfile.readline()
            tofile.write(seq)