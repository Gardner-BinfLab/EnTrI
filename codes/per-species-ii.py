#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from shutil import rmtree
from re import match

def makedir(dirname):
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot'
outdir = '/home/fatemeh/EnTrI/results/per-species-ii'
makedir(outdir)
species_names = ["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t"]
insertion_indices = []
list_of_files = listdir(clusters)
for filename in list_of_files:
    with open('{0}/{1}'.format(clusters, filename)) as from_file:
        for line in from_file:
            cells = line.split()
            insertion_indices.append((cells[1], cells[4]))
insertion_indices = list(set(insertion_indices))
insertion_indices.sort()

prev_name = ''
for item in insertion_indices:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item[0])
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item[0])
    if match_result:
        name = match_result.group(1)
    if name in species_names:
        if name != prev_name:
            if 'writefile' in locals():
                writefile.close()
            writefile = open(outdir + '/' + name + '.fasta', 'w')
        writefile.write(item[0] + '\t' + item[1] + '\n')
    prev_name = name