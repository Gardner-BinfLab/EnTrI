#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from os import listdir, mkdir, path
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


def add_insertion_index(in_clusters, in_iis, output):
    iis_dict = {}
    with open(in_iis) as iis_file:
        for line in iis_file:
            cells = line.split()
            iis_dict[cells[0]] = cells[1]
    list_of_files = listdir(in_clusters)
    for filename in list_of_files:
        with open('{0}/{1}'.format(in_clusters, filename)) as from_file:
            with open('{0}/{1}'.format(output, filename), 'w') as to_file:
                for line in from_file:
                    cells = line.split()
                    for item in cells:
                        to_file.write(item+'\t')
                    if cells[1] in iis_dict.keys():
                        to_file.write(iis_dict[cells[1]]+'\n')
                    else:
                        to_file.write('-1\n')

indir = str(argv[1])
iifile = str(argv[2])
outdir = str(argv[3])
makedir(outdir)
add_insertion_index(indir, iifile, outdir)