#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from argparse import ArgumentParser
from os import path, listdir, mkdir
from re import match, compile
from shutil import rmtree

def makedir(dirname):
    # print ('Making directory \'{0}\'..'.format(dirname))
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

parser = ArgumentParser(description='Amends names of genes in the given folder')
parser.add_argument('indir', help='Input directory (EFam-MSAs)')
parser.add_argument('clustiidir', help='Directory of clusters with coordinates and insertion indices (clust2plot_output)')
parser.add_argument('outdir', help='Directory for outputs')
args = parser.parse_args()
indir = args.indir
clustiidir = args.clustiidir
outdir = args.outdir
makedir(outdir)

gene_dict = {"BN373": "KlePn.2", "CS17": "EscCo.1", "Ecoli9000q": "EscCo.2", "ENC": "EntCo.1",
             "ERS227112": "KlePn.1", "ETEC": "EscCo.3", "NCTC13441": "EscCo.4", "ROD": "CitRo.1",
             "SEN": "SalEn.1", "SL1344": "SalEn.3", "STM": "SalEn.2", "STMMW": "SalEn.4",
             "t": "SalEn.5", "b": "EscCo.5"}

list_of_files = listdir(indir)
for filename in list_of_files:
    with open('{0}/{1}'.format(indir, filename)) as from_file:
        clustiiname = filename[0:-3]+'txt'
        with open('{0}/{1}'.format(clustiidir, clustiiname)) as clustii_file:
            with open('{0}/{1}'.format(outdir, filename), 'w') as to_file:
                for line in from_file:
                    if line.startswith('>'):
                        clustii_line = clustii_file.readline()
                        # print(clustii_line)
                        cells = clustii_line.split()
                        ii = float(cells[4])
                        if ii < 0:
                            essentiality =''
                        elif ii < 0.2:
                            essentiality = '_Ess'
                        elif ii < 2:
                            essentiality = '_NEs'
                        else:
                            essentiality = '_BeL'
                        #match_result = match('>\s*([A-Za-z]+(\d+_|_)?)(misc_RNA|R|_)?(\d+_\d+\-\d+)', line)
                        #line = line.replace(match_result.group(4), match_result.group(4)+essentiality)
                        #match_result = match('([a-zA-Z0-9]*)', match_result.group(1))
                        match_result = match('>\s*([a-zA-Z0-9]+?)_(\S*)', line)
                        pattern = compile('^\d+\-\d+\S*$')
                        if pattern.match(match_result.group(2)):
                            match_result = match('([a-zA-Z]*)\S*', match_result.group(1))
                        if match_result.group(1) in gene_dict.keys():
                            gene_name = gene_dict[match_result.group(1)]
                            line = line.replace(match_result.group(1), gene_name, 1)
                            match_result = match('>\s*\S*(\d+\-\d+)', line)
                            line = line.replace(match_result.group(1), match_result.group(1)+essentiality)
                    to_file.write(line)