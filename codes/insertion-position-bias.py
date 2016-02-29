from collections import defaultdict
from re import match
from os import listdir, path
from math import floor

seqdb = '../sequences/fasta-dna/chromosome/seqdb.fasta'
plots = '../sequences/plot-files/chromosome'
result = '../results/insertion-position-bias.out'
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
    genome_insertions[item] = sum(plots_dict[item])

gene_name = ''
position_insertion = [0]*180
with open(result, 'w') as tofile:
    with open(seqdb, 'r') as sequencefile:
        for line in sequencefile:
            if line.startswith('>'):
                if gene_name != '':
                    tofile.write(gene_name)
                    for item in position_insertion:
                        tofile.write('\t' + str(item))
                    tofile.write('\t' + str(ii) + '\n')
                match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is not None:
                    gene_name = match_result.group(1)
                    strain_name = match_result.group(2)
                    start = int(match_result.group(3)) - 1
                    end = int(match_result.group(4)) - 1
                    strand = match_result.group(5)
                else:
                    strain_name = ''
                if strain_name not in plots_dict.keys() or (end - start + 1) < 180:
                    gene_name = ''
                else:
                    for i in range(0,60):
                        if plots_dict[strain_name][start+i] > 0:
                            position_insertion[i] = float(1) / (float(genome_insertions[strain_name])/genome_length[strain_name])
                        else:
                            position_insertion[i] = float(0)
                        if plots_dict[strain_name][end-59+i] > 0:
                            position_insertion[120+i] = float(1) / (float(genome_insertions[strain_name])/genome_length[strain_name])
                        else:
                            position_insertion[120+i] = float(0)
                        interval = (end - 60) - (start + 60) + 1
                        ithstart = int(start + 60 + floor(float(interval * i)/60))
                        ithend = int(start + 60 + floor(float(interval * (i+1))/60))
                        length = ithend - ithstart
                        position_insertion[60+i] = float(float(sum([1 for x in plots_dict[strain_name][ithstart:ithend] if x > 0]))/length)/(float(genome_insertions[strain_name])/genome_length[strain_name])
                        gene_insertions = sum([1 for x in plots_dict[strain_name][start:end+1] if x > 0])
                        ii = (float(gene_insertions)/(end-start+1))/(float(genome_insertions[strain_name])/genome_length[strain_name])
                    if strand == 'Complement':
                        position_insertion.reverse()
