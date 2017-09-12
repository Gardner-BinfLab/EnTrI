from collections import defaultdict
from re import match
from os import listdir, path
from math import floor
from Bio import SeqIO


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict


seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
plots = '../data/plot-files/chromosome'
essentialitydir = '../results/biases/dbscan/'
outdir = '../results/5prime-starts.txt'

lfiles = listdir(plots)
lfiles.remove('BW25113.plot')
lfiles.remove('EC958.plot')
sequences = read_fasta_sequences(seqdb)

ins_st = 0
ins_nst = 0
nins_st = 0
nins_nst = 0
ins_trnst = 0


for item in lfiles:
    ins_dict = defaultdict(list)
    name, extension = path.splitext(item)
    with open('{0}/{1}'.format(plots, item), 'r') as plotfile:
        for line in plotfile:
            cells = line.split()
            summation = cells[0] + cells[1]
            ins_dict[name].append(int(int(cells[0]) > 0))

    list_of_files = listdir(essentialitydir)
    essentiality_dict = dict()
    for filename in list_of_files:
        with open(essentialitydir + filename) as fromfile:
            for line in fromfile:
                cells = line.split()
                essentiality_dict[cells[0]] = cells[2]

    with open(seqdb, 'r') as sequencefile:
        for line in sequencefile:
            if line.startswith('>'):
                match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is not None:
                    gene_name = match_result.group(1)
                    strain_name = match_result.group(2)
                    start = int(match_result.group(3)) - 1
                    end = int(match_result.group(4)) - 1
                    strand = match_result.group(5)

                    if strain_name in ins_dict.keys() and gene_name in essentiality_dict.keys() \
                            and essentiality_dict[gene_name] == 'essential':
                        if strand == 'Forward':
                            fiveperc = round((end - start + 1) * 5 / 100) + start
                            tenperc = round((end - start + 1) * 10 / 100)
                            if tenperc % 3:
                                tenperc += (3 - (tenperc % 3))
                            indices = [i for i, x in enumerate(ins_dict[strain_name][start:fiveperc + 1]) if x == 1]
                        elif strand == 'Complement':
                            fiveperc = end - round((end - start + 1) * 5 / 100)
                            tenperc = round((end - start + 1) * 10 / 100)
                            if tenperc % 3:
                                tenperc += (3 - (tenperc % 3))
                            indices = [i for i, x in enumerate(ins_dict[strain_name][fiveperc:end + 1]) if x == 1]
                            indices = [(end - fiveperc + 1) - i for i in indices]

                        frame2 = 0
                        if len(indices) > 0:
                            last_insertion = max(indices)
                            if last_insertion % 3 == 2:
                                frame2 = 1
                            if last_insertion % 3:
                                last_insertion += (3 - (last_insertion % 3))
                            flag = 0
                            for i in range(last_insertion, tenperc, 3):
                                if sequences[gene_name].seq[i:i + 3] == 'atg' or sequences[gene_name].seq[
                                                                                 i:i + 3] == 'ttg' or sequences[
                                                                                                          gene_name].seq[
                                                                                                      i:i + 3] == 'ctg' or \
                                                sequences[gene_name].seq[i:i + 3] == 'gtg':
                                    flag = 1
                            if flag:
                                ins_st += 1
                            else:
                                ins_nst += 1
                                if frame2 == 1:
                                    ins_trnst += 1
                        else:
                            flag = 0
                            for i in range(3, tenperc, 3):
                                if sequences[gene_name].seq[i:i + 3] == 'atg' or sequences[gene_name].seq[
                                                                                 i:i + 3] == 'ttg' or sequences[
                                                                                                          gene_name].seq[
                                                                                                      i:i + 3] == 'ctg' or \
                                                sequences[gene_name].seq[i:i + 3] == 'gtg':
                                    flag = 1
                            if flag:
                                nins_st += 1
                            else:
                                nins_nst += 1

                        # if len(indices) > 0:
                        #     last_insertion = max(indices)
                        #     if last_insertion % 3 == 2:
                        #         last_insertion += 1
                        #         flag = 0
                        #         for i in range(last_insertion, tenperc, 3):
                        #             if sequences[gene_name].seq[i:i + 3] == 'atg' or sequences[gene_name].seq[
                        #                                                          i:i + 3] == 'ttg' or sequences[
                        #                                                                                   gene_name].seq[
                        #                                                                               i:i + 3] == 'ctg' or \
                        #                         sequences[gene_name].seq[i:i + 3] == 'gtg':
                        #                 flag = 1
                        #         if flag:
                        #             ins_st += 1
                        #         else:
                        #             ins_nst += 1
                        # else:
                        #     flag = 0
                        #     for i in range(3, tenperc, 3):
                        #         if sequences[gene_name].seq[i:i + 3] == 'atg' or sequences[gene_name].seq[
                        #                                                          i:i + 3] == 'ttg' or sequences[
                        #                                                                                   gene_name].seq[
                        #                                                                               i:i + 3] == 'ctg' or \
                        #                         sequences[gene_name].seq[i:i + 3] == 'gtg':
                        #             flag = 1
                        #     if flag:
                        #         nins_st += 1
                        #     else:
                        #         nins_nst += 1

print(ins_st)
print(ins_nst)
print(nins_st)
print(nins_nst)
print(ins_trnst)