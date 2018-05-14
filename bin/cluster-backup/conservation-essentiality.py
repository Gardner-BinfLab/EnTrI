from os import listdir
conspath = '/data/people/fas31/homology/hmmsearch/'
# conspath = '/home/fatemeh/endosymbionts/hmmsearch/'
consdict = {}
list_of_files = listdir(conspath)
for filename in list_of_files:
    with open(conspath + filename, 'r') as fromfile:
        clusts = []
        for line in fromfile:
            cells = line.split()
            if len(cells) > 3 and cells[3].startswith('clust'):
                clusts.append(cells[3])
    clusts = list(set(clusts))
    for item in clusts:
        if item in consdict.keys():
            consdict[item] += 1
        else:
            consdict[item] = 1


writepath = '/data/people/fas31/homology/homologs.count'
# writepath = '/home/fatemeh/endosymbionts/homologs.count'
with open(writepath, 'w') as tofile:
    for key in consdict.keys():
        tofile.write(key+'\t'+str(consdict[key])+'\n')
    tofile.write('total'+'\t'+str(len(list_of_files))+'\n')
