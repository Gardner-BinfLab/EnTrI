from sys import argv
from os import path, mkdir, system, listdir, remove, getcwd
from shutil import rmtree, copyfile, move
from Bio import SeqIO
from collections import defaultdict
from re import match
from argparse import ArgumentParser


def makedir(dirname):
    # print ('Making directory \'{0}\'..'.format(dirname))
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Overwrite? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict


outdir = '../results/interesting_genes/'
makedir(outdir)
clustersfile = '../results/hieranoid/clusters.txt'
seqdbfile = '../data/fasta-protein/chromosome/seqdb.fasta'
essentialitydir = '../results/biases/dbscan'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
numstrains = 16
seqdb = read_fasta_sequences(seqdbfile)
essentiality = {i: 0 for i in seqdb.keys()}
npeq = {i: 4.5 for i in seqdb.keys()}
ii = {i: 7 for i in seqdb.keys()}

list_of_files = listdir(essentialitydir)
for filename in list_of_files:
    with open (essentialitydir + '/' + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            npeq[cells[0]] = cells[3]
            ii[cells[0]] = cells[1]
            if cells[2] == 'essential':
                essentiality[cells[0]] = 1

with open(k12path, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        essentiality[cells[0]] = 1
        npeq[cells[0]] = -4.5
        ii[cells[0]] = 0

species_names = {"salmonella": ["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"],
                 "ecoli": ["CS17", "ETEC", "NCTC13441", "b", "BW25113", "EC958"],
                 "klebsiella": ["ERS227112", "BN373"], "citrobacter": ["ROD"], "enterobacter": ["ENC"]}

interesting_genes = outdir + 'universally-conserved_always-essential.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t'+
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if len(cells) == numstrains and sum([essentiality[i] for i in cells]) == numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

interesting_genes = outdir + 'universally-conserved_sometimes-essential.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t' +
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if len(cells) == numstrains and 0 < sum([essentiality[i] for i in cells]) < numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

interesting_genes = outdir + 'universally-conserved.tsv'
start = '0'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t'+
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\tEsCoMG1655 position\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if len(cells) == numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                        if match_result.group(1) == 'b':
                            start = match('[a-zA-Z]+\S+\s+\[[^/]+/(\d+)', seqdb[item].description).group(1)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (ii[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\t' + start)
                tofile.write('\n')

interesting_genes = outdir + 'universally-conserved-log-ii-thresh.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t'+
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if len(cells) == numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

# interesting_genes = outdir + 'universally-conserved_universally-non-essential.tsv'
# with open(clustersfile, 'r') as fromfile:
#     with open(interesting_genes, 'w') as tofile:
#         tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
#                      'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
#                      'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
#                      'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
#                      'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
#                      'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
#                      'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t'+
#                      'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
#         for line in fromfile:
#             essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
#                                     "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
#                                     "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
#                                     "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
#                                     "ROD": ('X', 'X'), "ENC": ('X', 'X')}
#             cells = line.split()
#             if len(cells) == numstrains and sum([essentiality[i] for i in cells]) == 0:
#                 for item in cells:
#                     match_result = match('(\S+)_\S+', item)
#                     if not match_result:
#                         match_result = match('([a-zA-Z]+)\S+', item)
#                     strain = match_result.group(1)
#                     essentiality_strains[strain] = (npeq[item], item)
#                 product = ''
#                 gene = ''
#                 counter = 0
#                 while (product == '' or gene == '') and counter < len(cells):
#                     if product == '':
#                         product = seqdb[cells[counter]].description.split('] [')[1]
#                     if gene == '':
#                         gene = seqdb[cells[counter]].description.split('] [')[3]
#                     counter += 1
#                 tofile.write(gene + '\t' + product + '\t')
#                 tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
#                              str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
#                              '\t' +
#                              str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
#                              '\t' +
#                              str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
#                              str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
#                              str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
#                              '\t' +
#                              str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
#                              str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
#                 tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
#                              essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
#                              essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
#                              essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
#                              essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
#                              essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
#                              essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
#                              essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
#                 tofile.write('\n')

interesting_genes = outdir + 'universally-unconserved_always-essential.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t' +
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if sum([essentiality[i] for i in cells]) == len(cells) < numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

interesting_genes = outdir + 'universally-unconserved_sometimes-essential.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t' +
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if 0 < sum([essentiality[i] for i in cells]) < len(cells) < numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

interesting_genes = outdir + 'essential-in-some-species.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t' +
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()

            newcells = list(cells)
            for i in range(len(newcells)-1, -1, -1):
                if float(npeq[newcells[i]]) < 1.644854:
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
                    length.append(len(species_names[key]))

            if flag and len(newcells) == sum(length) and sum(length) < numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)

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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')

interesting_genes = outdir + 'sometimes-essential.tsv'
with open(clustersfile, 'r') as fromfile:
    with open(interesting_genes, 'w') as tofile:
        tofile.write('Gene\tProduct\tEnClNCTC9394 NPEQ\tKlPnEcl8 NPEQ\tKlPnRH201207 NPEQ\tCiRoICC168 NPEQ\t' +
                     'SaTySL1344 NPEQ\tSaTySL3261 NPEQ\tSaTyD23580 NPEQ\tSaTyA130 NPEQ\tSaEnP125109 NPEQ\t' +
                     'SaTyTy2 NPEQ\tEsCoEC958 NPEQ\tEsCoST131 NPEQ\tEsCoCS17 NPEQ\tEsCoH10407 NPEQ\t' +
                     'EsCoBW25113 NPEQ\tEsCoMG1655 NPEQ\tEnClNCTC9394 locus tag\tKlPnEcl8 locus tag\t' +
                     'KlPnRH201207 locus tag\tCiRoICC168 locus tag\tSaTySL1344 locus tag\tSaTySL3261 locus tag\t' +
                     'SaTyD23580 locus tag\tSaTyA130 locus tag\tSaEnP125109 locus tag\tSaTyTy2 locus tag\t' +
                     'EsCoEC958 locus tag\tEsCoST131 locus tag\tEsCoCS17 locus tag\tEsCoH10407 locus tag\t' +
                     'EsCoBW25113 locus tag\tEsCoMG1655 locus tag\n')
        for line in fromfile:
            essentiality_strains = {"SEN": ('X', 'X'), "SL1344": ('X', 'X'), "STM": ('X', 'X'), "STMMW": ('X', 'X'),
                                    "t": ('X', 'X'), "SL3261": ('X', 'X'), "CS17": ('X', 'X'), "ETEC": ('X', 'X'),
                                    "NCTC13441": ('X', 'X'), "b": ('X', 'X'), "BW25113": ('X', 'X'),
                                    "EC958": ('X', 'X'), "ERS227112": ('X', 'X'), "BN373": ('X', 'X'),
                                    "ROD": ('X', 'X'), "ENC": ('X', 'X')}
            cells = line.split()
            if 0 < sum([essentiality[i] for i in cells]) < numstrains:
                for item in cells:
                    match_result = match('(\S+)_\S+', item)
                    if not match_result:
                        match_result = match('([a-zA-Z]+)\S+', item)
                    strain = match_result.group(1)
                    essentiality_strains[strain] = (npeq[item], item)
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
                tofile.write(str(essentiality_strains['ENC'][0]) + '\t' + str(essentiality_strains['BN373'][0]) + '\t' +
                             str(essentiality_strains['ERS227112'][0]) + '\t' + str(essentiality_strains['ROD'][0]) +
                             '\t' +
                             str(essentiality_strains['SL1344'][0]) + '\t' + str(essentiality_strains['SL3261'][0]) +
                             '\t' +
                             str(essentiality_strains['STMMW'][0]) + '\t' + str(essentiality_strains['STM'][0]) + '\t' +
                             str(essentiality_strains['SEN'][0]) + '\t' + str(essentiality_strains['t'][0]) + '\t' +
                             str(essentiality_strains['EC958'][0]) + '\t' + str(essentiality_strains['NCTC13441'][0]) +
                             '\t' +
                             str(essentiality_strains['CS17'][0]) + '\t' + str(essentiality_strains['ETEC'][0]) + '\t' +
                             str(essentiality_strains['BW25113'][0]) + '\t' + str(essentiality_strains['b'][0]) + '\t')
                tofile.write(essentiality_strains['ENC'][1] + '\t' + essentiality_strains['BN373'][1] + '\t' +
                             essentiality_strains['ERS227112'][1] + '\t' + essentiality_strains['ROD'][1] + '\t' +
                             essentiality_strains['SL1344'][1] + '\t' + essentiality_strains['SL3261'][1] + '\t' +
                             essentiality_strains['STMMW'][1] + '\t' + essentiality_strains['STM'][1] + '\t' +
                             essentiality_strains['SEN'][1] + '\t' + essentiality_strains['t'][1] + '\t' +
                             essentiality_strains['EC958'][1] + '\t' + essentiality_strains['NCTC13441'][1] + '\t' +
                             essentiality_strains['CS17'][1] + '\t' + essentiality_strains['ETEC'][1] + '\t' +
                             essentiality_strains['BW25113'][1] + '\t' + essentiality_strains['b'][1])
                tofile.write('\n')