from collections import defaultdict, Counter
from re import match
from os import listdir, path
from math import floor
from Bio import SeqIO
from statistics import mode, median, mean
import matplotlib.pyplot as plt

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
plots = '../data/plot-files/chromosome'
core_essentials = '../results/define-core-accessory-hieranoid-fitch-core-essentials/'
outdir = '../results/5prime-ancestrally-essential.out'

list_of_files = listdir(plots)
plots_dict = defaultdict(list)
for filename in list_of_files:
    name, extension = path.splitext(filename)
    with open('{0}/{1}'.format(plots, filename), 'r') as plotfile:
        for line in plotfile:
            cells = line.split()
            plots_dict[name].append(int(cells[0]) + int(cells[1]))

list_of_files = listdir(core_essentials)
list_of_files.remove('alignments')
list_of_files.remove('profiles')
list_of_files.remove('hhmake')
list_of_files.remove('clusters.info')
essentiality_dict = {}
for filename in list_of_files:
    essentiality_dict[filename] = []
    with open(core_essentials + filename) as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                gene = line.rstrip()[1:]
                essentiality_dict[filename].append(gene)

sequences = read_fasta_sequences(seqdb)
clusts_dict = {key: [] for key in essentiality_dict}

for clust in essentiality_dict.keys():
    for gene in essentiality_dict[clust]:
        desc = sequences[gene].description
        match_result = match('\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', desc)
        if match_result is None:
            match_result = match('\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', desc)
        if match_result is not None:
            strain_name = match_result.group(2)
            start = int(match_result.group(3)) - 1
            end = int(match_result.group(4)) - 1
            strand = match_result.group(5)
            if strand == 'Forward':
                fiveperc = round((end - start + 1) * 5 / 100) + start + 1
                if sum(plots_dict[strain_name][start:fiveperc]) > 0:
                    clusts_dict[clust].append(gene)
            elif strand == 'Complement':
                fiveperc = round((end - start + 1) * 5 / 100) - end
                if sum(plots_dict[strain_name][fiveperc:end+1]) > 0:
                    clusts_dict[clust].append(gene)

counts = []
with open(outdir, 'w') as tofile:
    for key in clusts_dict.keys():
        tofile.write(key)
        counter = 0
        for item in clusts_dict[key]:
            tofile.write('\t' + item)
            counter += 1
        tofile.write('\t' + str(len(clusts_dict[key])) + '\n')
        counts.append(counter)

print(str(mean(counts)))
print(str(median(counts)))
print(str(mode(counts)))
hist = Counter(counts)
xvals = list(hist.keys())
xvals.sort()
yvals = [hist[item] for item in xvals]
plt.bar(xvals,height=yvals)
plt.show()