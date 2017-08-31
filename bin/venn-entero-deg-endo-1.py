from os import listdir, remove, system, makedirs
from shutil import copyfile
from itertools import combinations

address = {'deg': '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid-fitch/all/core-essential-genomes/hhmake/',
           'prot': '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid-fitch/proteobacteria/core-essential-genomes/hhmake/',
           'gamma': '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid-fitch/gammaproteobacteria/core-essential-genomes/hhmake/',
           'endo': '/home/fatemeh/EnTrI/results/endosymbionts/define-core-accessory-hieranoid/core-essential-genomes/hhmake/',
           'entero': '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials/hhmake/'}
outdir = '/home/fatemeh/EnTrI/results/venn-entero-deg-endo/'

list_of_degs = listdir(address['deg'])
list_of_prots = listdir(address['prot'])
list_of_gammas = listdir(address['gamma'])
list_of_endos = listdir(address['endo'])
list_of_enteros = listdir(address['entero'])

dict_of_lists = {'deg':list_of_degs, 'prot':list_of_prots, 'gamma':list_of_gammas, 'endo':list_of_endos, 'entero':list_of_enteros}
for key in dict_of_lists.keys():
    temppath = address[key] + 'temp.txt'
    for filename in dict_of_lists[key]:
        filepath = address[key] + filename
        copyfile(filepath, temppath)
        with open(temppath, 'r') as fromfile:
            with open(filepath, 'w') as tofile:
                for line in fromfile:
                    cells = line.split()
                    if len(cells) and cells[0] == 'NAME':
                        cells[1] = filename.split('.')[0]
                        line = '  '.join(cells) + '\n'
                    tofile.write(line)
    remove(temppath)
    system('cat ' + address[key] + '* > ' + address[key] + key + 'fam.hmm')

combs = list(combinations(dict_of_lists.keys(), 2))
# combs = [('deg','endo'), ('deg', 'gamma'), ('endo', 'gamma'), ('entero','deg'), ('entero', 'endo'), ('entero', 'gamma'),
#          ('prot','deg'), ('prot', 'endo'), ('prot','entero'),('prot','gamma')]

for subset in combs:
    combpath = outdir + subset[0] + subset[1] + '/'
    makedirs(combpath)
    for filename in dict_of_lists[subset[0]]:
        system('hhsearch -i ' + address[subset[0]] + filename + ' -d ' + address[subset[1]] + subset[1] + 'fam.hmm -o ' + combpath + filename)

dict_of_pairs = {}
for subset in combs:
    pair = subset[0] + subset[1]
    dict_of_pairs[pair] = {}
    combpath = outdir + pair + '/'
    for filename in listdir(combpath):
        with open(combpath + filename, 'r') as fromfile:
            firstclust = filename.split('.')[0]
            for line in fromfile:
                cells = line.split()
                if len(cells) > 1 and cells[1].startswith(subset[1] + '-clust'):
                    secondclust = cells[1]
                    prob = float(cells[2])
                    if prob == 100:
                        if secondclust not in dict_of_pairs[pair].keys():
                            dict_of_pairs[pair][secondclust] = (firstclust, prob)
                        elif prob > dict_of_pairs[pair][secondclust][1]:
                            dict_of_pairs[pair][secondclust] = (firstclust, prob)
                    break

clusts = []
for pair in dict_of_pairs.keys():
    for clust1 in dict_of_pairs[pair].keys():
        clust2 = dict_of_pairs[pair][clust1][0]
        clusts.append((clust1,clust2))

counter = 1
clust_nums = [0] * len(clusts)
for i in range(0, len(clusts)):
    if clust_nums[i] == 0:
        pair = clusts[i]
        clust_nums[i] = counter
        members = {pair[0], pair[1]}
        members_copy = {'empty'}
        while members != members_copy:
            members_copy = members
            for j in range(0, len(clusts)):
                if clusts[j][0] in members or clusts[j][1] in members:
                    members = members | {clusts[j][0], clusts[j][1]}
                    clust_nums[j] = counter
        counter += 1

clust_dict = {}
for i in range(0, len(clusts)):
    num = clust_nums[i]
    clust_dict['clust-' + str(num)] = [clusts[i][0], clusts[i][1]]
    for j in range(0, len(clusts)):
        if clust_nums[j] == num:
            clust_dict['clust-' + str(num)].append(clusts[j][0])
            clust_dict['clust-' + str(num)].append(clusts[j][1])


with open(outdir + 'final-clusters.txt', 'w') as tofile:
    for i in range(1, counter):
        key = 'clust-' + str(i)
        tofile.write(key)
        for item in list(set(clust_dict[key])):
            tofile.write('\t' + item)
        tofile.write('\n')


maxlen = 0
minlen = 6
for key in clust_dict.keys():
    length = len(list(set(clust_dict[key])))
    maxlen = max(maxlen, length)
    minlen = min(minlen, length)
    if length > 5:
        print(clust_dict[key])

print(maxlen)
print(minlen)