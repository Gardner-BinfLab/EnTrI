from collections import defaultdict
from re import match
from os import listdir, path
from math import floor

seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
plots = '../data/plot-files/chromosome'
result = '../results/insertion-indices/insertion-position-bias-first100.out'
essentialitydir = '../results/biases/dbscan/'
genome_length = {"SL1344":4878012, "STMMW":4879400, "SEN":4685848, "t":4791961, "STM":4895639, "ETEC":5153435,
                 "b":4641652, "CS17":4994793, "NCTC13441":5174631, "ROD":5346659, "BN373":5324709, "ERS227112":5869288,
                 "ENC":4908759, "SL3261":4878012, "EC958":5109767, "BW25113":4631469}
list_of_files = listdir(plots)
plots_dict = defaultdict(list)
for filename in list_of_files:
    name, extension = path.splitext(filename)
    with open('{0}/{1}'.format(plots, filename), 'r') as plotfile:
        for line in plotfile:
            cells = line.split()
            plots_dict[name].append(int(cells[0]) + int(cells[1]))

list_of_files = listdir(essentialitydir)
essentiality_dict = dict()
for filename in list_of_files:
    with open(essentialitydir + filename) as fromfile:
        for line in fromfile:
            cells = line.split()
            essentiality_dict[cells[0]] = cells[2]

genome_insertions = dict()
for item in plots_dict.keys():
    genome_insertions[item] = sum([1 for x in plots_dict[item] if x > 0])

position_insertion = [0] * 200
sums = [0]*200
with open(result, 'w') as tofile:
    with open(seqdb, 'r') as sequencefile:
        for line in sequencefile:
            if line.startswith('>'):
                match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is not None:
                    gene_name = match_result.group(1)
                    strain_name = match_result.group(2)
                    start = int(match_result.group(3)) - 1
                    end = int(match_result.group(4)) - 1
                    strand = match_result.group(5)
                    if strain_name in plots_dict.keys() and (end - start + 1) >= 200 and\
                                    gene_name in essentiality_dict.keys() and essentiality_dict[gene_name]=='essential':
                        for i in range(0, 200):
                            position_insertion[i] = int(plots_dict[strain_name][start+i] > 0)
                        if strand == 'Complement':
                            temp = position_insertion[0:(len(position_insertion)-1)]
                            temp.reverse()
                            temp.append(position_insertion[len(position_insertion)-1])
                            position_insertion = temp
                        tofile.write(gene_name)
                        for item in position_insertion:
                            tofile.write('\t' + str(item))
                        tofile.write('\n')
                        sums = [x + y for x, y in zip(sums, position_insertion)]
print(sums)