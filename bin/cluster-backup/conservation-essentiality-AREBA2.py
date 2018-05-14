from os import listdir
from collections import defaultdict
from re import match

distfile = '/data/people/fas31/homology/dnadist.out'
names = dict()
distances = defaultdict(list)
counter = 0
strain = ''
with open(distfile, 'r') as fromfile:
    line = fromfile.readline()
    for line in fromfile:
        cells = line.split()
        if match('[a-zA-Z]+\d+',cells[0]):
            strain = cells[0]
            names[strain] = counter
            counter += 1
            distances[strain] = cells[1:]
        else:
            distances[strain].extend(cells)

clusts = defaultdict(list)
hmmerpath = '/data/people/fas31/homology/hmmsearch/'
list_of_files = listdir(hmmerpath)
for filename in list_of_files:
    strain = filename.split('.')[0]
    with open(hmmerpath + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            if len(cells) > 3 and cells[3].startswith('clust'):
                clusts[cells[3]].append(strain)

for key in clusts.keys():
    clusts[key] = list(set(clusts[key]))

writepath = '/data/people/fas31/homology/conservation2.out'
with open(writepath, 'w') as tofile:
    for key in clusts.keys():
        conservation= 0
        for indx, item1 in enumerate(clusts[key][:-1]):
            for item2 in clusts[key][indx+1:]:
                # if item1 in distances.keys() and item2 in distances.keys():
                try:
                    num = names[item2]
                    dist = float(distances[item1][num])
                    conservation = max(conservation, dist)
                except:
                    pass
        # print(key+'\t'+str(conservation))
        tofile.write(key+'\t'+str(conservation)+'\n')

