#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from re import match, findall
from os import listdir, mkdir, path, remove
from Bio import SeqIO
from shutil import rmtree
from collections import defaultdict


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


def add_coordinates(sequences, input_dir, output_dir):
    fastas = read_fasta_sequences(sequences)
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as from_file:
            with open('{0}/{1}'.format(output_dir, filename), 'w') as to_file:
                for line in from_file:
                    start = []
                    end = []
                    prot_start = []
                    prot_end = []
                    cells = line.split()
                    match_result = match('\S+\s+\[\S+/([\d\-\s]+)\((\S+)\)\]', fastas[cells[1]].description)
                    # if match_result is None: #Used for Matt's fasta file
                        # match_result = match('\S+\s+\S+\s+[^/]+/([\d\-\s]+):[\+\-]', fastas[cells[1]].description) #Used for Matt's fasta file
                    # strand = match_result.group(2)
                    find_result = findall('\d+\-\d+', match_result.group(1))
                    for item in find_result:
                        start_end = match('(\d+)\-(\d+)', item)
                        prot_start.append(int(start_end.group(1)))
                        prot_end.append(int(start_end.group(2)))
                    dom_start = int(cells[2])
                    dom_end = int(cells[3])
                    # check for sequences with join
                    # summation = 0
                    # i = 0
                    # while len(prot_end) > i+1 and (dom_start - 1) * 3 >=\
                    #         summation + prot_end[i+1] - prot_start[i+1] + 1:
                    #     i += 1
                    #     summation += prot_end[i] - prot_start[i] + 1
                    # start.append(prot_start[i] + (dom_start - 1) * 3 - summation)
                    # while len(prot_end) > i+1 and\
                    #         dom_end * 3 - 1 >= summation + prot_end[i+1] - prot_start[i+1] + 1:
                    #     end.append(prot_end[i])
                    #     i += 1
                    #     summation += prot_end[i] - prot_start[i] + 1
                    #     start.append(prot_start[i])
                    # end.append(prot_start[i] + dom_end * 3 - 1 - summation)
                    # write_start_end = ''
                    # for i in range(0, len(start)):
                    #     write_start_end += str(start[i])+'-'+str(end[i])
                    #     if i < len(start)-1:
                    #         write_start_end += ','
                    # write_to_file = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(cells[0], cells[1], cells[2], cells[3],
                    #                                                    write_start_end)
                    # to_file.write(write_to_file)
                    i = 0
                    while (dom_start - 1) * 3 >= prot_end[i] - prot_start[i] + 1:
                        dom_start -= int((prot_end[i] - prot_start[i] + 1) / 3)
                        dom_end -= int((prot_end[i] - prot_start[i] + 1) / 3)
                        i += 1
                    start.append(prot_start[i] + (dom_start - 1) * 3)
                    while len(prot_end) >= i+1 and dom_end * 3 - 1 >= prot_end[i] - prot_start[i] + 1:
                        end.append(prot_end[i])
                        dom_end -= int((prot_end[i] - prot_start[i] + 1) / 3)
                        i += 1
                        if len(prot_end) >= i+1:
                            start.append(prot_start[i])
                    if len(prot_end) >= i+1:
                        end.append(prot_start[i] + dom_end * 3 - 1)
                    write_start_end = ''
                    for i in range(0, len(start)):
                        write_start_end += str(start[i])+'-'+str(end[i])
                        if i < len(start)-1:
                            write_start_end += ','
                    write_to_file = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(cells[0], cells[1], cells[2], cells[3],
                                                                       write_start_end)
                    to_file.write(write_to_file)

def remove_empty_files(input_dir):
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        if not path.getsize('{0}/{1}'.format(input_dir, filename)):
            remove('{0}/{1}'.format(input_dir, filename))

def add_insertion_index(in_clusters, in_plots, output):
    list_of_files = listdir(in_plots)
    plots_dict = defaultdict(list)
    for filename in list_of_files:
        name, extension = path.splitext(filename)
        with open('{0}/{1}'.format(in_plots, filename)) as from_file:
            for line in from_file:
                cells = line.split()
                plots_dict[name].append(int(cells[0]) + int(cells[1]))
    list_of_files = listdir(in_clusters)
    for filename in list_of_files:
        with open('{0}/{1}'.format(in_clusters, filename)) as from_file:
            with open('{0}/{1}'.format(output, filename), 'w') as to_file:
                for line in from_file:
                    cells = line.split()
                    reads = 0
                    gene_length = 0
                    find_result = findall('\d+\-\d+', cells[4])
                    match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
                    if not match_result:
                        match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
                    if match_result:
                        name = match_result.group(1)
                        if name in plots_dict.keys():
                            genome_insertions = sum([1 for x in plots_dict[name] if x > 0])
                            genome_length = len(plots_dict[name])
                            for item in find_result:
                                start_end = match('(\d+)\-(\d+)', item)
                                start = int(start_end.group(1))
                                end = int(start_end.group(2))
                                reads += sum(plots_dict[name][start-1:end])
                                gene_insertions = sum([1 for x in plots_dict[name][start-1:end] if x > 0])
                                gene_length += end - start + 1
                            insertion_index = (float(gene_insertions) / gene_length) / (float(genome_insertions) /genome_length )
                            # insertion_index = float(gene_insertions) / gene_length # used for Matt
                        else:
                            reads = -1
                            insertion_index = -1
                            for item in find_result:
                                start_end = match('(\d+)\-(\d+)', item)
                                start = int(start_end.group(1))
                                end = int(start_end.group(2))
                                gene_length += end - start + 1
                        to_file.write('{0}\t{1}\t{2}\t{3}\n'.format(line.rstrip('\n'), reads, insertion_index, gene_length))
    remove_empty_files(output)

# seqdb = 'EnTrI-on-real-data/sequences/fasta/seqdb.fasta'
# indir = 'EnTrI-on-real-data/homclust-output/EFam-clusters'
# plotsdir = 'EnTrI-on-real-data/plot-files'
# outdir = 'clust2plottest'

seqdb = str(argv[1])
indir = str(argv[2])
plotsdir = str(argv[3])
outdir = str(argv[4])
makedir(outdir)
clusters = '{0}/clusters_with_coordinates'.format(outdir)
makedir(clusters)
add_coordinates(seqdb, indir, clusters)
plots = '{0}/final_clusters'.format(outdir)
makedir(plots)
add_insertion_index(clusters, plotsdir, plots)