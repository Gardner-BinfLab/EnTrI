from os import listdir
from Bio import SeqIO

cmsearch = '/data/people/fas31/homology/cmsearch/'
genomes = '/data/genomes/bacteria2016-04-14/fasta/'
outdir = '/data/people/fas31/homology/cmsearch-sequences.fa'

list_of_files = listdir(genomes)
list_of_cmsearch = listdir(cmsearch)
for filename in list_of_files:
    newfilename = filename.split('.')[0] + '.tblout'
    if newfilename not in list_of_cmsearch:
        list_of_files.remove(filename)
with open(outdir, 'w') as tofile:
    for filename in list_of_files:
        newfilename = filename.split('.')[0] + '.tblout'
        with open(genomes + filename, 'rU') as fasta_file:
            sequence = list(SeqIO.parse(fasta_file, 'fasta'))[0]
        with open(cmsearch+newfilename, 'r') as fromfile:
            for line in fromfile:
                if not line.startswith('#'):
                    cells = line.split()
                    start = int(cells[7])
                    end = int(cells[8])
                    strand = cells[9]
                    if strand == '+':
                        seq = sequence.seq[start-1:end]
                    else:
                        seq2 = sequence.seq[end-1:start][::-1]
                        seq = ''
                        for item in seq2:
                            if item in ['A','a']:
                                seq += 't'
                            elif item in ['T', 't', 'U', 'u']:
                                seq += 'a'
                            elif item in ['C', 'c']:
                                seq += 'g'
                            elif item in ['G', 'g']:
                                seq += 'c'
                            else:
                                seq += 'n'
                    tofile.write('>' + sequence.id[:-1] + '\n')
                    tofile.write(str(seq) + '\n')
                    break
