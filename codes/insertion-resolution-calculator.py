from collections import defaultdict
from re import match
from os import listdir, path

plots = '../sequences/plot-files/chromosome'
genome_length = {"SL1344":4878012, "STMMW":4879400, "SEN":4685848, "t":4791961, "STM":4895639, "ETEC":5153435,
                 "b":4641652, "CS17":4994793, "NCTC13441":5174631, "ROD":5346659, "BN373":5324709, "ERS227112":5869288,
                 "ENC":4908759}
list_of_files = listdir(plots)
plots_dict = defaultdict(list)
for filename in list_of_files:
    name, extension = path.splitext(filename)
    with open('{0}/{1}'.format(plots, filename), 'r') as plotfile:
        for line in plotfile:
            cells = line.split()
            plots_dict[name].append(int(cells[0]) + int(cells[1]))

genome_insertions = dict()
for item in plots_dict.keys():
    genome_insertions[item] = sum([1 for x in plots_dict[item] if x > 0])

for item in plots_dict.keys():
    print ('Genome name: ' + item + '\t' + 'Resolution: ' + str(float(genome_length[item])/float(genome_insertions[item])) + '\n')