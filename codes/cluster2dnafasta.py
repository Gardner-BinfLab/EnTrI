#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from os import path, mkdir, system, listdir, remove, getcwd
from shutil import rmtree, copyfile, move
from Bio import SeqIO
from collections import defaultdict
from re import search
from argparse import ArgumentParser

def check_directory_name(dirname):
    if dirname.endswith('/'):
        dirname = dirname[:-1]
    return dirname

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

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

parser = ArgumentParser(description='Get DNA fasta files of the clusters made by homclust')
#parser.add_argument('dupthr', help='Only clusters with more number of members than this threshold are considered')
parser.add_argument('seqdb', help='Fasta file for DNA sequences in all species')
parser.add_argument('cldir', help='Directory of EFam clusters resulted from homclust')
#parser.add_argument('msadir', help='Directory of EFam MSAs resulted from homclust')
parser.add_argument('outdir', help='Directory for outputs')
args = parser.parse_args()
#dup_thresh = int(args.dupthr)
seqdb = args.seqdb
cldir = args.cldir
outdir = args.outdir
#msadir = args.msadir
cldir = check_directory_name(cldir)
#msa_indir = check_directory_name(msadir)
outdir = check_directory_name(outdir)
makedir(outdir)
fasta_outdir = outdir
#fasta_outdir = '{}/fasta'.format(outdir)
#makedir(fasta_outdir)
#msa_outdir = '{}/msa'.format(outdir)
#makedir('{}/msa'.format(outdir))

fastas = read_fasta_sequences(seqdb)
list_of_files = listdir(cldir)
for filename in list_of_files:
    num_lines = sum(1 for line in open('{0}/{1}'.format(cldir, filename)))
    #if num_lines >= dup_thresh:
    name, exten = path.splitext(filename)
    with open('{0}/{1}'.format(cldir, filename), 'r') as cluster_file:
        #with open('{0}/{1}.msa'.format(msa_indir, name), 'r') as msa_infile:
        with open('{0}/{1}.fasta'.format(fasta_outdir, name), 'w') as fasta_file:
            for line in cluster_file:
                cells = line.split()
                id = cells[1]
                efam_id = cells[0]
                start = (int(cells[2]) - 1) * 3
                end = int(cells[3]) * 3
                #fastas[name].id = fastas[id].id+'_{}-{}'.format(cells[2], cells[3])
                fasta_file.write('>')
                fasta_file.write(fastas[id].description)
                #search_result = search('EFam\-\d+_(\d+_\d+)', efam_id)
                #fasta_file.write(search_result.group(1))
                for i in range(start, end):
                    if not ((i - start) % 60):
                        fasta_file.write('\n')
                        #print (fastas[id].id, i, start, end, cells[2], cells[3])
                    fasta_file.write(fastas[id].seq[i])
                fasta_file.write('\n')

        # with open('{0}/{1}'.format(cldir, filename), 'r') as cluster_file:
        #     with open('{0}/{1}.msa'.format(msa_indir, name), 'r') as msa_infile:
        #         with open('{0}/{1}.msa'.format(msa_outdir, name), 'w') as msa_outfile:
        #             for line in msa_infile:
        #                 if line.startswith('>'):
        #                     efam_id = cluster_file.readline().split()[0]
        #                     msa_outfile.write('>')
        #                     search_result = search('EFam\-\d+_(\d+_\d+)', efam_id)
        #                     msa_outfile.write(search_result.group(1))
        #                     msa_outfile.write('\n')
        #                 else:
        #                     msa_outfile.write(line)