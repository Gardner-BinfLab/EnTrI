from re import findall
from math import floor, ceil
from Bio import SeqIO
from itertools import product

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
                    if (sum(list(map(float, cells[1:21]))) / 20) > (sum(list(map(float, cells[21:81]))) / 60) + 0.02:
                        tofile.write('>' + cells[0])
                        tofile.write('\n')
                        fiveperc = ceil(floor(len(str(sequence_dict[cells[0]].seq))*20 / 100)/3)*3
                        tofile.write(str(sequence_dict[cells[0]].seq)[0:fiveperc])
                        tofile.write('\n')
                    elif (sum(list(map(float, cells[1:21]))) / 20) < (sum(list(map(float, cells[21:81]))) / 60) - 0.02:
                        tofilenins.write('>' + cells[0])
                        tofilenins.write('\n')
                        fiveperc = ceil(floor(len(str(sequence_dict[cells[0]].seq)) * 20 / 100)/3)*3
                        tofilenins.write(str(sequence_dict[cells[0]].seq)[0:fiveperc])
                        tofilenins.write('\n')

list_of_codons = [''.join(p) for p in product('actg', repeat=3)]
worddictmore = dict((el,0) for el in list_of_codons)
worddictless = dict((el,0) for el in list_of_codons)

codoncountermore = 0
with open(resultins, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('>'):
            codons = findall('[atgc]{3}', line[3:len(line)])
            codoncountermore += len(codons)
            codons = list(set(codons))
            for item in codons:
                worddictmore[item] += 1

codoncounterless = 0
with open(resultnins, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('>'):
            codons = findall('[atgc]{3}', line[3:len(line)])
            codoncounterless += len(codons)
            codons = list(set(codons))
            for item in codons:
                worddictless[item] += 1

with open(info, 'w') as tofile:
    tofile.write('codon\tMoreThanAve\tLessThanAve\n')
    for item in list_of_codons:
        tofile.write(item+'\t'+str(worddictmore[item])+'\t'+str(worddictless[item])+'\n')
    tofile.write('sum\t'+str(codoncountermore)+'\t'+str(codoncounterless))
    # tofile.write('sum\t' + str(sum(worddictmore.values())) + '\t' + str(sum(worddictless.values())))