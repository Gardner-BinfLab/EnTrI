from os import listdir
from Bio import SeqIO
from re import match

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

clustersfile = '../results/hieranoid/clusters.txt'
seqdbfile = '../data/fasta-protein/chromosome/seqdb.fasta'
essentialitydir = '../results/biases/normalised-pca'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
intereseting_genes = '../results/interesting_genes.csv'

seqdb = read_fasta_sequences(seqdbfile)
essentiality = {i: 0 for i in seqdb.keys()}
ii = {i: 0 for i in seqdb.keys()}

list_of_files = listdir(essentialitydir)
for filename in list_of_files:
    with open (essentialitydir + '/' + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            ii[cells[0]] = cells[1]
            if cells[2] == 'essential':
                essentiality[cells[0]] = 1

with open(k12path, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        essentiality[cells[0]] = 1
        ii[cells[0]] = 7

species_names = {"salmonella":["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"],
    "ecoli":["CS17", "ETEC", "NCTC13441", "b","BW25113", "EC958"], "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "enterobacter":
    ["ENC"]}

with open(clustersfile, 'r') as fromfile:
    with open(intereseting_genes, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            # if (len(cells) >= 2 and 0 < sum([essentiality[i] for i in cells]) < len(cells)) or (len(cells) < 16 and
            # 0 < sum([essentiality[i] for i in cells])): # all 500 genes with differentiating essentiality

            # if (len(cells) == 16 and 14 < sum([essentiality[i] for i in cells]) < 16) or\
            #         (14 < len(cells) < 16 and 14 < sum([essentiality[i] for i in cells])): # all clusters with only 1 exception

            newcells = list(cells)
            for i in range(len(newcells)-1, -1, -1):
                if float(ii[newcells[i]]) < 1.644854:
                    del newcells[i]
            for i in range(len(newcells)):
                match_result = match('(\S+)_\S+', newcells[i])
                if not match_result:
                    match_result = match('([a-zA-Z]+)\S+', newcells[i])
                newcells[i] = match_result.group(1)

            flag = 0
            length = []
            for key in species_names.keys():
                if set(species_names[key]) <= set(newcells):
                    flag = 1
                    # length = max(length, len(species_names[key]))
                    length.append(len(species_names[key]))

            if flag and len(newcells) == sum(length) and sum(length) < 16:
            # if flag and len(newcells) <= 9:
                product = ''
                gene = ''
                counter = 0
                while (product == '' or gene == '') and counter < len(cells):
                    if product == '':
                        product = seqdb[cells[counter]].description.split('] [')[1]
                    if gene == '':
                        gene = seqdb[cells[counter]].description.split('] [')[3]
                    counter += 1
                tofile.write(gene + '\t' + product + '\t')
                for item in cells:
                    tofile.write(item + '\t' + str(ii[item]) + '\t')
                tofile.write('\n')
