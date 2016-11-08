from os import listdir
from collections import defaultdict
import numpy as np

fastas_dir = '../data/fasta-protein/chromosome'
plots_dir = '../data/plot-files/chromosome'
list_of_files = listdir(plots_dir)
plots = defaultdict(list)
sampledplots = []
sumlength = []
weights = []
numsamples = 1000
locusid = []
list_of_files = ['t.plot']

for filename in list_of_files:
    plotfile = []
    with open(plots_dir + '/' + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            plotfile.append(cells)
    lid = filename.split('.')[0]
    locusid.append(lid)
    plots[lid] = [int(plotfile[i][0])+int(plotfile[i][1]) for i in range(len(plotfile))]
    for i in range(numsamples):
        sampledplots.append(np.random.choice(plots[lid], len(plots[lid])))