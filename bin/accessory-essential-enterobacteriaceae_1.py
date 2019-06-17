from os import listdir
ecoli = []
ecolidir = '../results/define-core-accessory-hieranoid-fitch/ecoli/core-essential-genomes/'
ecolis = listdir(ecolidir)
for filename in ecolis:
    with open(ecolidir + filename, 'r') as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                cells = line.split()
                ecoli.append(cells[0][1:])

citrobacter = []
citrobacterdir = '../results/define-core-accessory-hieranoid-fitch/citrobacter/core-essential-genomes/'
citrobacters = listdir(citrobacterdir)
for filename in citrobacters:
    with open(citrobacterdir + filename, 'r') as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                cells = line.split()
                citrobacter.append(cells[0][1:])

klebsiella = []
klebsielladir = '../results/define-core-accessory-hieranoid-fitch/klebsiella/core-essential-genomes/'
klebsiellas = listdir(klebsielladir)
for filename in klebsiellas:
    with open(klebsielladir + filename, 'r') as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                cells = line.split()
                klebsiella.append(cells[0][1:])

salmonella = []
salmonelladir = '../results/define-core-accessory-hieranoid-fitch/salmonella/core-essential-genomes/'
salmonellas = listdir(salmonelladir)
for filename in salmonellas:
    with open(salmonelladir + filename, 'r') as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                cells = line.split()
                salmonella.append(cells[0][1:])

# enterobacter = []
# enterobacterdir = '../results/define-core-accessory-hieranoid-fitch/enterobacter/core-essential-genomes/'
# enterobacters = listdir(enterobacterdir)
# for filename in enterobacters:
#     with open(enterobacterdir + filename, 'r') as fromfile:
#         for line in fromfile:
#             if line.startswith('>'):
#                 cells = line.split()
#                 enterobacter.append(cells[0][1:])
#
# all = ecoli + citrobacter + salmonella + enterobacter + klebsiella
all = ecoli + citrobacter + salmonella + klebsiella

eggnog = '../results/eggnog-mapper/seqdb.fasta.emapper.annotations'
eggnogdict = {}
cogdict = {}
with open(eggnog, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        if cells[4] != '' and not cells[0].startswith('ENC'):
            eggnogdict[cells[0]] = cells[4]
        elif not cells[0].startswith('ENC'):
            whole = cells[8].split(',')
            for item in whole:
                pre = item.split('@')[0]
                post = item.split('@')[1]
                if post == 'NOG':
                    cogdict[cells[0]] = pre

genedict = {}
clusters = '../results/hieranoid/clusters.txt'
with open(clusters, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        cells = [item for item in cells if not item.startswith('ENC')]
        if len(set(cells) & set(all)) > 0:
            flag = 0
            for item in set(cells) & set(all):
                if item in eggnogdict.keys():
                    name = eggnogdict[item]
                    flag = 1
                    break
            if not flag:
                for item in set(cells) & set(all):
                    if item in eggnogdict.keys():
                        name = cogdict[item]
                        flag = 1
                        break
            if not flag:
                print(item)
            else:
                genedict[name] = [0, 0, 0, 0]
            if len(set(cells) & set(citrobacter)) > 0:
                genedict[name][0] = 1
            if len(set(cells) & set(salmonella)) > 0:
                genedict[name][1] = 1
            if len(set(cells) & set(ecoli)) > 0:
                genedict[name][2] = 1
            if len(set(cells) & set(klebsiella)) > 0:
                genedict[name][3] = 1
            # if len(set(cells) & set(enterobacter)) > 0:
            #     genedict[name][4] = 1

for item in genedict.keys():
    genedict[item] = tuple(genedict[item])

sorted_dict = sorted(genedict, key=genedict.get, reverse=True)
writepath = '../results/define-core-accessory-hieranoid-fitch/heatmap.tsv'
with open(writepath, 'w') as tofile:
    tofile.write('gene\tCitrobacter\tSalmonella\tEscherichia\tKlebsiella\n')
    for item in sorted_dict:
        first = item[0:3]
        last = item[3:]
        gene = first.lower() + last
        tofile.write(gene)
        for i in genedict[item]:
            tofile.write('\t' + str(i))
        tofile.write('\n')