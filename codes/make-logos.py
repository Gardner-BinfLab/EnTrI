#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, remove

genomespath = '../sequences/fasta-genome/chromosome'
plotspath = '../sequences/plot-files/chromosome'
genomextension = 'fa'
logospath = '../results/logos.txt'

if path.isfile(logospath):
    input_var = 'i'
    while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
        input_var = input('File {0} already exists. Replace? [Y/n] '.format(logospath))
        if input_var == 'Y' or input_var == 'y':
            remove(logospath)
        else:
            raise SystemExit

list_of_files = listdir(plotspath)
for filename in list_of_files:
    plots = []
    sequences = ''
    name, extension = path.splitext(filename)
    with open('{0}/{1}'.format(plotspath, filename)) as plotfile:
        for line in plotfile:
            cells = line.split()
            plots.append(int(cells[0]) + int(cells[1]))
    with open('{0}/{1}.{2}'.format(genomespath, name, genomextension)) as genomefile:
        for line in genomefile:
            if not line.startswith('>'):
                if line.endswith('\n'):
                    line = line[:-1]
                sequences += line
    with open (logospath, 'a') as logosfile:
        for i in range(0, len(plots)):
            if plots[i] > 0:
                logosfile.write('>{0}_{1}\n'.format(name, i+1))
                if 10 <= i < len(sequences)-10:
                    logosfile.write('{0}\n'.format(sequences[i-10:i+11]))
                elif i < 10:
                    logosfile.write('{0}\n'.format(sequences[(i-10)%len(sequences):len(sequences)]+sequences[0:i+11]))
                else:
                    logosfile.write('{0}\n'.format(sequences[i-10:len(sequences)]+sequences[0:(i+11)%len(sequences)]))