from collections import defaultdict
from re import match
from os import listdir, path

seqdb = '../sequences/fasta-dna/chromosome/seqdb.fasta'
plots = '../sequences/plot-files/chromosome'
result = '../results/check-biases.out'
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

gene_name = ''
with open(result, 'w') as tofile:
    with open(seqdb, 'r') as sequencefile:
        for line in sequencefile:
            if line.startswith('>'):
                if gene_name != '':
                    tofile.write('{0}\t{1}\t{2}\t{3}\n'.format(gene_name, ii, float(start)/genome_length[strain_name], float(gc)/(end-start+1)))
                gc = 0
                match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\(', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\(', line)
                if match_result is not None:
                    gene_name = match_result.group(1)
                    strain_name = match_result.group(2)
                    start = int(match_result.group(3))
                    end = int(match_result.group(4))
                else:
                    strain_name = ''
                if strain_name not in plots_dict.keys():
                    gene_name = ''
                else:
                    gene_insertions = sum([1 for x in plots_dict[strain_name][start-1:end] if x > 0])
                    ii = (float(gene_insertions)/(end-start+1))/(float(genome_insertions[strain_name])/genome_length[strain_name])
            else:
                gc += line.count('g') + line.count('c')
