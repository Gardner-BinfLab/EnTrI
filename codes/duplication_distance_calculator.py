#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from Bio import Phylo
from os import listdir, system, path, mkdir, remove, stat
from shutil import rmtree
from argparse import ArgumentParser
from re import match
from csv import reader

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

parser = ArgumentParser(description='Plots distance vs. delta insertion index')
parser.add_argument('msadir', help='Directory of nucleotide multiple sequence alignments (pal2nal)')
parser.add_argument('iidir', help='Directory of clusters with insertion indices (clust2plot_output/final_clusters)')
parser.add_argument('distdir', help='Directory for saving the distances between sequences')
parser.add_argument('deltaiidir', help='Directory for saving the distance between insertion indices')
args = parser.parse_args()
msadir = args.msadir
iidir = args.iidir
distdir = args.distdir
deltaiidir = args.deltaiidir

makedir(distdir)
makedir(deltaiidir)

list_of_files = listdir(msadir)
for filename in list_of_files:
    if stat('{0}/{1}'.format(msadir, filename)).st_size:
        try:
            system('fdnadist -sequence {0}/{1} -method F84 -outfile {2}/{3}.temp'.format(msadir, filename, distdir, path.splitext(filename)[0]))
        except:
            raise SystemExit
        with open('{0}/{1}'.format(msadir, filename)) as msafile:
            with open('{0}/{1}.temp'.format(distdir, path.splitext(filename)[0])) as distfile_from:
                with open('{0}/{1}.txt'.format(distdir, path.splitext(filename)[0]), 'w') as distfile_to:
                    flag = 1
                    distfile_from.readline()
                    while True:
                        msaline = msafile.readline()
                        while msaline and not msaline.startswith('>'):
                            msaline = msafile.readline()
                        if not msaline:
                            break
                        match_result = match('>(\S*)', msaline)
                        replacement = match_result.group(1)
                        line_to_write = distfile_from.readline()
                        line_to_write = line_to_write.rstrip()
                        match_result = match('([a-zA-Z])+\S*', line_to_write)
                        while not match_result:
                            cells = line_to_write.split()
                            line_to_write = '\t'+'\t'.join(cells)
                            distfile_to.write(line_to_write)
                            line_to_write = distfile_from.readline()
                            line_to_write = line_to_write.rstrip()
                            match_result = match('([a-zA-Z])+\S*', line_to_write)
                        cells = line_to_write.split()
                        cells[0] = replacement
                        if not flag:
                            line_to_write = '\n'+'\t'.join(cells)
                        else:
                            flag = 0
                            line_to_write = '\t'.join(cells)
                        distfile_to.write(line_to_write)
                    line_to_write = distfile_from.readline()
                    while line_to_write:
                        line_to_write = line_to_write.rstrip()
                        cells = line_to_write.split()
                        line_to_write = '\t'+'\t'.join(cells)
                        distfile_to.write(line_to_write)
                        line_to_write = distfile_from.readline()
        remove('{0}/{1}.temp'.format(distdir, path.splitext(filename)[0]))

list_of_files = listdir(iidir)
for filename in list_of_files:
    with open('{0}/{1}'.format(iidir, filename)) as iifile:
        rows = reader(iifile, delimiter='\t')
        iitable = []
        for row in rows:
            iitable.append(row)
    flag = 0
    with open('{0}/{1}.txt'.format(deltaiidir, path.splitext(filename)[0]), 'w') as deltaiifile:
        for row1 in iitable:
            if flag:
                deltaiifile.write('\n')
            flag = 1
            deltaiifile.write('{0}_{1}-{2}'.format(row1[1],row1[2],row1[3]))
            for row2 in iitable:
                if float(row1[6]) == -1 or float(row2[6])== -1:
                    deltaii = -1
                elif max(float(row1[6]), float(row2[6])) > 0:
                    deltaii = abs(float(row1[6]) - float(row2[6]))
                    # deltaii = abs(float(row1[6]) - float(row2[6])) / max(float(row1[6]), float(row2[6]))
                else:
                    deltaii = 0
                deltaiifile.write('\t{0}'.format(deltaii))
        deltaiifile.write('\n')