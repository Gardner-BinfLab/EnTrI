from Bio import SeqIO
from re import match, findall
from collections import defaultdict
from os import listdir


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

rootpath = '../results/define-core-accessory-hieranoid/all/core-essential-genomes'
ecolipath = '../results/define-core-accessory-hieranoid/ecoli/core-essential-genomes'
salmonellapath = '../results/define-core-accessory-hieranoid/salmonella/core-essential-genomes'
klebsiellapath = '../results/define-core-accessory-hieranoid/klebsiella/core-essential-genomes'
klebsiellaenterobacterpath = '../results/define-core-accessory-hieranoid/klebsiellaenterobacter/core-essential-genomes'
salmonellacitrobacterpath = '../results/define-core-accessory-hieranoid/salmonellacitrobacter/core-essential-genomes'
salmonellaecolicitrobacterpath = '../results/define-core-accessory-hieranoid/salmonellaecolicitrobacter/core-essential-genomes'
genomepath = '../data/fasta-protein/chomosome/seqdb.fasta'

ecodes = '../results/branch-root/root-ecoli-description.txt'
necodes = '../results/branch-root/not-root-ecoli-description.txt'
ecoinfo = '../results/branch-root/ecoli-info.txt'

saldes = '../results/branch-root/root-salmonella-description.txt'
nsaldes = '../results/branch-root/not-root-salmonella-description.txt'
salinfo = '../results/branch-root/salmonella-info.txt'

kledes = '../results/branch-root/root-klebsiella-description.txt'
nkledes = '../results/branch-root/not-root-klebsiella-description.txt'
kleinfo = '../results/branch-root/klebsiella-info.txt'

kleentdes = '../results/branch-root/root-klebsiellaenterobacter-description.txt'
nkleentdes = '../results/branch-root/not-root-klebsiellaenterobacter-description.txt'
kleentinfo = '../results/branch-root/klebsiellaenterobacter-info.txt'

salcitdes = '../results/branch-root/root-salmonellacitrobacter-description.txt'
nsalcitdes = '../results/branch-root/not-root-salmonellacitrobacter-description.txt'
salcitinfo = '../results/branch-root/salmonellacitrobacter-info.txt'

salecocitdes = '../results/branch-root/root-salmonellaecolicitrobacter-description.txt'
nsalecocitdes = '../results/branch-root/not-root-salmonellaecolicitrobacter-description.txt'
salecocitinfo = '../results/branch-root/salmonellaecolicitrobacter-info.txt'

list_of_files = listdir(rootpath)
alldict = {}
for filename in list_of_files:
    alldict.update(read_fasta_sequences(rootpath+'/'+filename))

list_of_files = listdir(ecolipath)
ecodict = {}
for filename in list_of_files:
    ecodict.update(read_fasta_sequences(ecolipath+'/'+filename))

list_of_files = listdir(salmonellapath)
saldict = {}
for filename in list_of_files:
    saldict.update(read_fasta_sequences(salmonellapath+'/'+filename))

list_of_files = listdir(klebsiellapath)
kledict = {}
for filename in list_of_files:
    kledict.update(read_fasta_sequences(klebsiellapath+'/'+filename))

list_of_files = listdir(klebsiellaenterobacterpath)
kleentdict = {}
for filename in list_of_files:
    kleentdict.update(read_fasta_sequences(klebsiellaenterobacterpath+'/'+filename))

list_of_files = listdir(salmonellacitrobacterpath)
salcitdict = {}
for filename in list_of_files:
    salcitdict.update(read_fasta_sequences(salmonellacitrobacterpath+'/'+filename))

list_of_files = listdir(salmonellaecolicitrobacterpath)
salecocitdict = {}
for filename in list_of_files:
    salecocitdict.update(read_fasta_sequences(salmonellaecolicitrobacterpath+'/'+filename))

with open(necodes, 'w') as nfile:
    with open(ecodes, 'w') as yfile:
        for item in sorted(ecodict.keys()):
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?", ecodict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item+'\t'+match_result.group(1)+'\t'+match_result.group(2)+'\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(ecodes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(necodes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(ecoinfo, 'w') as tofile:
    tofile.write('WORD\tEcoli\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(nsaldes, 'w') as nfile:
    with open(saldes, 'w') as yfile:
        for item in sorted(saldict.keys()):
            match_result = match("\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?", saldict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(saldes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nsaldes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(salinfo, 'w') as tofile:
    tofile.write('WORD\tSalmonella\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(nkledes, 'w') as nfile:
    with open(kledes, 'w') as yfile:
        for item in sorted(kledict.keys()):
            match_result = match(
                "\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?",
                kledict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(kledes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nkledes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(kleinfo, 'w') as tofile:
    tofile.write('WORD\tKlebsiella\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(nkleentdes, 'w') as nfile:
    with open(kleentdes, 'w') as yfile:
        for item in sorted(kleentdict.keys()):
            match_result = match(
                "\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?",
                kleentdict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(kleentdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nkleentdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(kleentinfo, 'w') as tofile:
    tofile.write('WORD\tKlebsiellaenterobacter\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(nsalcitdes, 'w') as nfile:
    with open(salcitdes, 'w') as yfile:
        for item in sorted(salcitdict.keys()):
            match_result = match(
                "\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?",
                salcitdict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(salcitdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nsalcitdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(salcitinfo, 'w') as tofile:
    tofile.write('WORD\tSalmonellacitrobacter\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))


with open(nsalecocitdes, 'w') as nfile:
    with open(salecocitdes, 'w') as yfile:
        for item in sorted(salecocitdict.keys()):
            match_result = match(
                "\S+\s+\[\S+\s+\([a-zA-Z]+\)\]\s+\[([^\]]+)\]\s+(?:\[[^\]]+\])?\s+(?:\[([a-zA-Z0-9]+)\])?",
                salecocitdict[item].description)
            if match_result:
                if item in alldict.keys():
                    if match_result.group(2):
                        yfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        yfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')
                else:
                    if match_result.group(2):
                        nfile.write(item + '\t' + match_result.group(1) + '\t' + match_result.group(2) + '\n')
                    else:
                        nfile.write(item + '\t' + match_result.group(1) + '\tNoname' + '\n')

worddict = defaultdict(list)

with open(salecocitdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][0] += 1
                else:
                    worddict[word] = [1, 0]

with open(nsalecocitdes, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        description = cells[1].replace(',', ' ')
        description = description.replace(':', ' ')
        description = description.replace('_', ' ')
        description = description.replace(';', ' ')
        #description = description.replace('-', ' ')
        description = description.replace('(', ' ')
        description = description.replace(')', ' ')
        description = description.replace('.', ' ')
        description = description.replace('/', ' ')
        description = description.replace('%', ' ')
        description = description.replace('*', ' ')
        description = description.replace('\\', ' ')
        description = description.replace('\'', ' ')
        description = description.replace('[', ' ')
        description = description.replace(']', ' ')
        words = findall("\S+", description)
        for item in words:
            if len(item) > 3:
                word = item.lower()
                if word in worddict.keys():
                    worddict[word][1] += 1
                else:
                    worddict[word] = [0, 1]

with open(salecocitinfo, 'w') as tofile:
    tofile.write('WORD\tSalmonellaecolicitrobacter\tNot\n')
    for item in sorted(worddict, key=worddict.get, reverse=True):
        tofile.write(item + '\t' + str(worddict[item][0]) + '\t' + str(worddict[item][1]) + '\n')
    sum_of_lists = [sum(i) for i in zip(*worddict.values())]
    tofile.write('SUM\t' + str(sum_of_lists[0]) + '\t' + str(sum_of_lists[1]))