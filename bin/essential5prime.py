from collections import defaultdict
from re import findall
from os import listdir, path
from math import floor, ceil
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
positional_insertions = '../results/insertion-indices/insertion-position-bias.out'
resultins = '../results/insertion-indices/essential5prime-insertions.out'
resultnins = '../results/insertion-indices/essential5prime-noinsertions.out'
info = '../results/insertion-indices/essential5prime-info.out'

sequence_dict = read_fasta_sequences(seqdb)

with open(positional_insertions, 'r') as fromfile:
    with open(resultins, 'w') as tofile:
        with open(resultnins, 'w') as tofilenins:
            for line in fromfile:
                cells = line.split()
                if float(cells[101]) < 0.2:
                    if (sum(list(map(float, cells[1:6]))) / 5) > (sum(list(map(float, cells[6:81]))) / 75) + 0.02:
                        tofile.write('>' + cells[0])
                        tofile.write('\n')
                        fiveperc = ceil(floor(len(str(sequence_dict[cells[0]].seq))*5 / 100)/3)*3
                        tofile.write(str(sequence_dict[cells[0]].seq)[0:fiveperc])
                        tofile.write('\n')
                    elif (sum(list(map(float, cells[1:6]))) / 5) < (sum(list(map(float, cells[6:81]))) / 75) - 0.02:
                        tofilenins.write('>' + cells[0])
                        tofilenins.write('\n')
                        fiveperc = ceil(floor(len(str(sequence_dict[cells[0]].seq)) * 5 / 100)/3)*3
                        tofilenins.write(str(sequence_dict[cells[0]].seq)[0:fiveperc])
                        tofilenins.write('\n')


with open(info, 'w') as tofile:
    worddict = {}
    with open(resultins, 'r') as fromfile:
        for line in fromfile:
            if not line.startswith('>'):
                codons = findall('[atgcATGC]{3}', line[3:len(line)])
                for item in codons:
                    codon = item.lower()
                    if codon in worddict.keys():
                        worddict[codon] += 1
                    else:
                        worddict[codon] = 1
    tofile.write('##################### More insertions than average #####################\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item]) + '\n')
    tofile.write('Sum\t'+ str(sum(worddict.values())) + '\n')

    worddict = {}
    with open(resultnins, 'r') as fromfile:
        for line in fromfile:
            if not line.startswith('>'):
                codons = findall('[atgcATGC]{3}', line[3:len(line)])
                for item in codons:
                    codon = item.lower()
                    if codon in worddict.keys():
                        worddict[codon] += 1
                    else:
                        worddict[codon] = 1
    tofile.write('##################### Less insertions than average #####################\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item]) + '\n')
    tofile.write('Sum\t' + str(sum(worddict.values())) + '\n')