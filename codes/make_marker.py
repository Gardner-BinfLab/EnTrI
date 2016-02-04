#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from argparse import ArgumentParser
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

parser = ArgumentParser(description='Finds marker genes to make phylogenetic trees')
parser.add_argument('indir', help='Directory of fasta clusters (EFam-fastas)')
parser.add_argument('indnadir', help='Directory of DNA fasta clusters (dna-clusters)')
parser.add_argument('outdir', help='Directory for outputs (protein fastas)')
parser.add_argument('outdnadir', help='Directory for outputs(DNA fastas)')
args = parser.parse_args()
input_dir = args.indir
input_dnadir = args.indnadir
output_dir = args.outdir
output_dnadir = args.outdnadir

makedir(output_dir)
makedir(output_dnadir)

gene_dict = {"BN373": 0, "CS17": 0, "ENC": 0, "ERS227112": 0, "ETEC": 0, "NCTC13441": 0,
             "ROD": 0, "SEN": 0, "SL1344": 0, "STM": 0, "STMMW": 0, "t": 0, "b":0}

name_dict = {"BN373": "KlePn.2", "CS17": "EscCo.1", "ENC": "EntCo.1",
             "ERS227112": "KlePn.1", "ETEC": "EscCo.2", "NCTC13441": "EscCo.3", "ROD": "CitRo.1",
             "SEN": "SalEn.1", "SL1344": "SalEn.3", "STM": "SalEn.2", "STMMW": "SalEn.4",
             "t": "SalEn.5", "b": "EscCo.4"}

list_of_gene_files = []

list_of_files = listdir(input_dir)
for filename in list_of_files:
    for key in gene_dict.keys():
        gene_dict[key] = 0
    with open('{0}/{1}'.format(input_dir, filename)) as from_file:
        for line in from_file:
            if line.startswith('>'):
                match_result = match('>\s*([a-zA-Z0-9]+?)_(\S*)', line)
                pattern = compile('^\d+\-\d+\S*$')
                if pattern.match(match_result.group(2)):
                    match_result = match('([a-zA-Z]*)\S*', match_result.group(1))
                if match_result.group(1) in gene_dict.keys():
                    gene_dict[match_result.group(1)] += 1
    if gene_dict[min(gene_dict, key=gene_dict.get)] == 1 and gene_dict[max(gene_dict, key=gene_dict.get)] == 1:
        list_of_gene_files.append(path.basename(filename))
#print (list_of_gene_files)

for key in name_dict.keys():
    open('{0}/{1}.fasta'.format(output_dir,name_dict[key]), 'a').close()
    open('{0}/{1}.fasta'.format(output_dnadir,name_dict[key]), 'a').close()

for filename in list_of_gene_files:
    with open('{0}/{1}'.format(input_dir, filename)) as from_file:
        line = from_file.readline()
        while line:
            if line.startswith('>'):
                match_result = match('>\s*([a-zA-Z0-9]+?)_(\S*)', line)
                pattern = compile('^\d+\-\d+\S*$')
                if pattern.match(match_result.group(2)):
                    match_result = match('([a-zA-Z]*)\S*', match_result.group(1))
                if match_result.group(1) in name_dict.keys():
                    with open('{0}/{1}.fasta'.format(output_dir, name_dict[match_result.group(1)]), 'a') as to_file:
                        to_file.write(line)
                        line = from_file.readline()
                        while line and not line.startswith('>'):
                            to_file.write(line)
                            line = from_file.readline()
            else:
                print("error")
                break

for filename in list_of_gene_files:
    with open('{0}/{1}'.format(input_dnadir, filename)) as from_file:
        line = from_file.readline()
        while line:
            if line.startswith('>'):
                match_result = match('>\s*([a-zA-Z0-9]+?)_\S+', line)
                if not match_result:
                    match_result = match('>\s*([a-zA-Z]+?)\d+', line)
                if match_result.group(1) in name_dict.keys():
                    with open('{0}/{1}.fasta'.format(output_dnadir, name_dict[match_result.group(1)]), 'a') as to_file:
                        to_file.write(line)
                        line = from_file.readline()
                        while line and not line.startswith('>'):
                            to_file.write(line)
                            line = from_file.readline()
            else:
                print("error")
                break

list_of_files = listdir(output_dir)
for filename in list_of_files:
    concatenated = ''
    with open('{0}/{1}'.format(output_dir, filename)) as from_file:
        for line in from_file:
            if not line.startswith('>'):
                concatenated += line.rstrip()
    with open('{0}/{1}'.format(output_dir, filename), 'w') as to_file:
        to_file.write('>{0}'.format(path.splitext(filename)[0]))
        counter = 0
        for item in concatenated:
            if not (counter % 60):
                to_file.write('\n')
            counter += 1
            to_file.write(item)

list_of_files = listdir(output_dnadir)
for filename in list_of_files:
    concatenated = ''
    with open('{0}/{1}'.format(output_dnadir, filename)) as from_file:
        for line in from_file:
            if not line.startswith('>'):
                concatenated += line.rstrip()
    with open('{0}/{1}'.format(output_dnadir, filename), 'w') as to_file:
        to_file.write('>{0}'.format(path.splitext(filename)[0]))
        counter = 0
        for item in concatenated:
            if not (counter % 60):
                to_file.write('\n')
            counter += 1
            to_file.write(item)
