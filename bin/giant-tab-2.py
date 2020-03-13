from os import listdir, path, mkdir
from shutil import rmtree
from re import match
import pandas as pd

def read_k12(inpath):
    all = []
    with open(inpath) as from_file:
        for line in from_file:
            cell = line.split()
            all.append(cell[0])
    return all

outdir = '../results/giant-tab/'
outpath = outdir + 'giant-tab-all.tsv'


tags = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008', 'DEG1006',
        'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002', 'DEG1001']

tradis = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'STM', 'STMMW', 'SL3261',
                'SL1344', 'SEN']
enterobacter = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
                'SL3261', 'SL1344', 'SEN']
endosymbiont = ['A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
                'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON',
                'BA000003', 'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG',
                'Sant', 'SOPEG']
gammaprot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003']
prot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
gammaprotexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b',
                   'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003']
protexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113',
              'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035',
              'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
bactexsymb = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029',
              'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
              'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020',
              'DEG1046', 'DEG1008', 'DEG1006', 'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002',
              'DEG1001']

lenenterobacter = len(enterobacter)
lenendosymbiont = len(endosymbiont)
lengammaprot = len(gammaprot)
lenprot = len(prot)
lengammaprotexsymb = len(gammaprotexsymb)
lenprotexsymb = len(protexsymb)
lenbact = len(tags)
lenbactexsymb = len(bactexsymb)

essentiality = '../results/biases/dbscan/'
ii_dict = dict()
list_of_files = listdir(essentiality)
for filename in list_of_files:
    with open(essentiality + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            ii_dict[cells[0]] = cells[3]

k12path = '../results/ecogene-k12.txt'
k12genes = read_k12(k12path)

clusters = '../results/all-hieranoid/clusters.txt'
with open(outpath, 'w') as tofile:
    tofile.write('Cluster\t')
    for item in tags[:-1]:
        tofile.write(item + '\t')
    tofile.write(tags[-1] + '\t')
    tofile.write('Enterobacteriaceae_80\tEnterobacteriaceae_90\tEnterobacteriaceae_100\t' +
                 'Endosymbiont_80\tEndosymbiont_90\tEndosymbiont_100\t' +
                 'Gammaproteobacter_80\tGammaproteobacter_90\tGammaproteobacter_100\t' +
                 'Gammaproteobacter-no-symb_80\tGammaproteobacter-no-symb_90\tGammaproteobacter-no-symb_100\t'
                 'Proteobacter_80\tProteobacter_90\tProteobacter_100\t' +
                 'Proteobacter-no-symb_80\tProteobacter-no-symb_90\tProteobacter-no-symb_100\t' +
                 'Bacteria_80\tBacteria_90\tBacteria_100\t' +
                 'Bacteria-no-symb_80\tBacteria-no-symb_90\tBacteria-no-symb_100\n')

    with open(clusters, 'r') as fromfile:
        cnt = 0
        for line in fromfile:
            genes = line.split()
            genes = [item for item in genes if not item.startswith('exDEG')]
            if len(genes) > 0:
                tofile.write('Clust' + str(cnt))
                clustdict = {el:0 for el in tags}
                cnt += 1

                numenterobacter = 0
                numendosymbiont = 0
                numgammaprot = 0
                numprot = 0
                numgammaprotexsymb = 0
                numprotexsymb = 0
                numbact = 0
                numbactexsymb = 0

                for gene in genes:
                    if gene.startswith('DEG'):
                        species = gene[0:7]
                    else:
                        species = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', str(gene)).group(1).strip('_')

                    if species in tradis:
                        if species != 'b':
                            if gene in ii_dict.keys():
                                clustdict[species] = int(float(ii_dict[gene]) <= 0)
                        else:
                            clustdict[species] = int(gene in k12genes)
                    else:
                        clustdict[species] = 1

                    if species in enterobacter and clustdict[species]:
                        numenterobacter += 1
                    if species in endosymbiont and clustdict[species]:
                        numendosymbiont += 1
                    if species in gammaprot and clustdict[species]:
                        numgammaprot += 1
                    if species in gammaprotexsymb and clustdict[species]:
                        numgammaprotexsymb += 1
                    if species in prot and clustdict[species]:
                        numprot += 1
                    if species in protexsymb and clustdict[species]:
                        numprotexsymb += 1
                    if species in tags and clustdict[species]:
                        numbact += 1
                    if species in bactexsymb and clustdict[species]:
                        numbactexsymb += 1

                for item in tags:
                    tofile.write('\t' + str(clustdict[item]))
                tofile.write('\t')
                tofile.write(str(int(numenterobacter * 100 / lenenterobacter >= 80)) + '\t')
                tofile.write(str(int(numenterobacter * 100 / lenenterobacter >= 90)) + '\t')
                tofile.write(str(int(numenterobacter * 100 / lenenterobacter >= 100)) + '\t')
                tofile.write(str(int(numendosymbiont * 100 / lenendosymbiont >= 80)) + '\t')
                tofile.write(str(int(numendosymbiont * 100 / lenendosymbiont >= 90)) + '\t')
                tofile.write(str(int(numendosymbiont * 100 / lenendosymbiont >= 100)) + '\t')
                tofile.write(str(int(numgammaprot * 100 / lengammaprot >= 80)) + '\t')
                tofile.write(str(int(numgammaprot * 100 / lengammaprot >= 90)) + '\t')
                tofile.write(str(int(numgammaprot * 100 / lengammaprot >= 100)) + '\t')
                tofile.write(str(int(numgammaprotexsymb * 100 / lengammaprotexsymb >= 80)) + '\t')
                tofile.write(str(int(numgammaprotexsymb * 100 / lengammaprotexsymb >= 90)) + '\t')
                tofile.write(str(int(numgammaprotexsymb * 100 / lengammaprotexsymb >= 100)) + '\t')
                tofile.write(str(int(numprot * 100 / lenprot >= 80)) + '\t')
                tofile.write(str(int(numprot * 100 / lenprot >= 90)) + '\t')
                tofile.write(str(int(numprot * 100 / lenprot >= 100)) + '\t')
                tofile.write(str(int(numprotexsymb * 100 / lenprotexsymb >= 80)) + '\t')
                tofile.write(str(int(numprotexsymb * 100 / lenprotexsymb >= 90)) + '\t')
                tofile.write(str(int(numprotexsymb * 100 / lenprotexsymb >= 100)) + '\t')
                tofile.write(str(int(numbact * 100 / lenbact >= 80)) + '\t')
                tofile.write(str(int(numbact * 100 / lenbact >= 90)) + '\t')
                tofile.write(str(int(numbact * 100 / lenbact >= 100)) + '\t')
                tofile.write(str(int(numbactexsymb * 100 / lenbactexsymb >= 80)) + '\t')
                tofile.write(str(int(numbactexsymb * 100 / lenbactexsymb >= 90)) + '\t')
                tofile.write(str(int(numbactexsymb * 100 / lenbactexsymb >= 100)) + '\n')

