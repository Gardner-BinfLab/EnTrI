#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from argparse import ArgumentParser
from os import path, listdir, mkdir
from re import match
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

gene_dict = {"KP1": "KlePn.1", "BN373": "KlePn.2", "NCTC13441": "EscCo.1", "CSSP291": "CroSa.1", "ENC": "EntCl.1",
             "ROD": "CitRo.1", "SO": "SheOn.1", "Sde": "SacDe.1", "PSE": "PseFO.1", "Gura": "GeoUr.1",
             "HRM2": "DesAu.1", "ELI": "EubLi.1", "BBR47": "BreBr.1", "BL": "BacLi.1", "Spirs": "SpiSm.1",
             "MM": "MetMa.1", "Natoc": "NatOc.1", "Dfer": "DyaFe.1", "SGRA": "SapGr.1", "t": "SalEn.1"}

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
                        ii = float(cells[6])
                        if ii < 0:
                            essentiality =''
                        elif ii < 0.2:
                            essentiality = '_Ess'
                        elif ii < 2:
                            essentiality = '_NEs'
                        else:
                            essentiality = '_BeL'
                        match_result = match('>\s*([A-Za-z]+(\d+_|_)?)(misc_RNA|R|_)?(\d+_\d+\-\d+)', line)
                        line = line.replace(match_result.group(4), match_result.group(4)+essentiality)
                        match_result = match('([a-zA-Z0-9]*)', match_result.group(1))
                        gene_name = gene_dict[match_result.group(1)]
                        line = line.replace(match_result.group(1), gene_name, 1)
                    to_file.write(line)