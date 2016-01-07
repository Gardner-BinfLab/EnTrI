from collections import defaultdict
from re import match
from os import listdir, path

seqdb = '../sequences/fasta-dna/chromosome/seqdb.fasta'
plots = '../sequences/plot-files/chromosome'
info = '../sequences/dnaa/dnaa-coordinates.txt'
result = '../results/check-biases.out'

info_table = defaultdict(list)
with open(info, 'r') as fastafile:
    for line in fastafile:
        cells = line.split()
        match_result = match('(\d+)\-(\d+)', cells[2])
        info_table[cells[0]]= [int(cells[1]), int(match_result.group(1)), int(match_result.group(2))]

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
                    tofile.write('{0}\t{1}\t{2}\t{3}\n'.format(gene_name, ii, float(distance)/info_table[strain_name][0], float(gc)/(end-start+1)))
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
                    ii = (float(gene_insertions)/(end-start+1))/(float(genome_insertions[strain_name])/info_table[strain_name][0])
                    dist1 = max(start, info_table[strain_name][1]) - min(end, info_table[strain_name][2]) - 1
                    dist2 = info_table[strain_name][0] - max(end, info_table[strain_name][2]) + min(start, info_table[strain_name][1]) -1
                    distance = min(dist1, dist2)
            else:
                gc += line.count('g') + line.count('c')
