#! /usr/bin/python3

##! /bin/sh
#""":"
#exec python3 $0 ${1+"$@"}
#"""
#__doc__ = """...Whatever..."""

from os import listdir, path, mkdir, system
from argparse import ArgumentParser
from shutil import rmtree, copy
from re import findall

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

parser = ArgumentParser(description='Reconciles gene tree and species tree')
parser.add_argument('msadir', help='Directory of MSA clusters (EFam-MSAs-edited-gene-names)')
parser.add_argument('dnadir', help='Directory of DNA fasta clusters (dna-clusters)')
parser.add_argument('speciestree', help='Species tree file (speciestree.fproml)')
parser.add_argument('pal2naldir', help='Directory of pal2nal outputs')
parser.add_argument('genetreedir', help='Directory of gene trees')
parser.add_argument('notungdir', help='Directory of Notung outputs')
args = parser.parse_args()
msadir = args.msadir
dnadir = args.dnadir
speciestree = args.speciestree
pal2naldir = args.pal2naldir
genetreedir = args.genetreedir
notungdir = args.notungdir
reconcile = '{0}/reconcile'.format(notungdir)
root =  '{0}/root'.format(notungdir)
png = '{0}/png'.format(notungdir)
png_to_study = '{0}/png_to_study'.format(notungdir)

makedir(pal2naldir)
makedir(genetreedir)
makedir(notungdir)
makedir(reconcile)
makedir(root)
makedir(png)
makedir(png_to_study)

list_of_files = listdir(dnadir)
for filename in list_of_files:
    try:
        system('pal2nal.pl {0}/{1}.msa {2}/{1}.fasta -output fasta > {3}/{1}.pal2nal'.format(msadir, path.splitext(filename)[0], dnadir, pal2naldir))
    except:
        raise SystemExit

list_of_files = listdir(pal2naldir)
for filename in list_of_files:
    try:
        system('FastTree -nosupport -quiet -nt {0}/{1}.pal2nal > {2}/{1}.newick'.format(pal2naldir, path.splitext(filename)[0], genetreedir))
    except:
        raise SystemExit

list_of_files = listdir(genetreedir)
for filename in list_of_files:
    try:
        system('java -jar ~/program-bank/Notung-2.6/Notung-2.6.jar -g {0}/{1} -s {2} --reconcile --nolosses --speciestag prefix --outputdir {3}/ --silent'.format(genetreedir, filename, speciestree, reconcile))
    except:
        raise SystemExit

list_of_files = listdir(reconcile)
for filename in list_of_files:
    try:
        system('java -jar ~/program-bank/Notung-2.6/Notung-2.6.jar {0}/{1} --root --nolosses --speciestag prefix --outputdir {2}/ --silent'.format(reconcile, filename, root))
    except:
        raise SystemExit

list_of_files = listdir(root)
for filename in list_of_files:
    try:
        system('java -jar ~/program-bank/Notung-2.6/Notung-2.6.jar {0}/{1} --savepng --nolosses --speciestag prefix --outputdir {2}/ --silent'.format(root, filename, png))
    except:
        raise SystemExit

list_of_files = listdir(root)
for filename in list_of_files:
    matches = []
    with open('{0}/{1}'.format(root, filename)) as from_file:
        for line in from_file:
            matches = matches + findall('Ess|NEs|BeL',line)
    ess = 0
    ness = 0
    for item in matches:
        if item == 'Ess':
            ess += 1
        else:
            ness += 1
    if ess > 2 and ness > 2:
        try:
            system('java -jar ~/program-bank/Notung-2.6/Notung-2.6.jar {0}/{1} --savepng --nolosses --speciestag prefix --outputdir {2}/ --silent'.format(root, filename, png_to_study))
        except:
            raise SystemExit