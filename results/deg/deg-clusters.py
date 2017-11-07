from re import search, match
annotations = 'DEG.annotations'
outpathdel = 'eggnog.clusters.del'
outpath = 'eggnog.clusters'
clusters = {}
with open(annotations, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        gene = cells[0]
        search_result = search(',(\S+)@bactNOG', cells[8])
        if not search_result:
            search_result = match('(\S+)@bactNOG', cells[8])
        if search_result:
            clust = search_result.group(1)
            if clust not in clusters.keys():
                clusters[clust] = [gene]
            else:
                clusters[clust].append(gene)

with open(outpathdel, 'w') as tofile:
    for key in clusters.keys():
        for item in clusters[key]:
            tofile.write(item + '\t')
        tofile.write('\n')

with open(outpathdel, 'r') as fromfile:
    with open(outpath, 'w') as tofile:
        for line in fromfile:
            genomes = {'DEG1001': [], 'DEG1003': [], 'DEG1006': [], 'DEG1008': [], 'DEG1013': [], 'DEG1015': [], 'DEG1020': [],
                       'DEG1023': [], 'DEG1027': [], 'DEG1029': [], 'DEG1034': [], 'DEG1036': [], 'DEG1038': [], 'DEG1040': [],
                       'DEG1042': [], 'DEG1045': [], 'DEG1047': [], 'DEG1002': [], 'DEG1005': [], 'DEG1007': [], 'DEG1012': [],
                       'DEG1014': [], 'DEG1017': [], 'DEG1021': [], 'DEG1024': [], 'DEG1028': [], 'DEG1031': [], 'DEG1035': [],
                       'DEG1037': [], 'DEG1039': [], 'DEG1041': [], 'DEG1043': [], 'DEG1046': []}
            cells = line.split()
            for item in cells:
                genomes[item[0:7]].append(item)
            while sum(len(genomes[i]) for i in genomes) > 0:
                for key in genomes.keys():
                    if len(genomes[key]) > 0:
                        item = genomes[key].pop()
                        tofile.write(item + '\t')
                tofile.write('\n')
