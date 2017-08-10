from re import match
dict = {}
with open('summary.txt', 'r') as fromfile:
    line = fromfile.readline()
    for line in fromfile:
        cells = line.split('\t')
        dict[cells[0]] = cells[5].replace(' ', '').lower()
level1 = 'bacteria;proteobacteria;gammaproteobacteria;enterobacteriales;enterobacteriaceae;'
level2 = 'bacteria;proteobacteria;gammaproteobacteria;enterobacteriales;'
level3 = 'bacteria;proteobacteria;gammaproteobacteria;'
level4 = 'bacteria;proteobacteria;'
level5 = 'bacteria;'

bacteria = []
proteobacteria = []
gammaproteobacteria = []

for key in dict.keys():
    if dict[key].startswith(level3) and not dict[key].startswith(level2):
        name = match(level3+'([^;]*)', dict[key]).group(1)
        mtc = 0
        for item in gammaproteobacteria:
            if dict[item].startswith(level3+name):
                mtc += 1
        if not mtc:
            gammaproteobacteria.append(key)

    elif dict[key].startswith(level4) and not dict[key].startswith(level3):
        name = match(level4+'([^;]*)', dict[key]).group(1)
        mtc = 0
        for item in proteobacteria:
            if dict[item].startswith(level4+name):
                mtc += 1
        if not mtc:
            proteobacteria.append(key)

    elif dict[key].startswith(level5) and not dict[key].startswith(level4):
        name = match(level5+'([^;]*)', dict[key]).group(1)
        mtc = 0
        for item in bacteria:
            if dict[item].startswith(level5+name):
                mtc += 1
        if not mtc:
            bacteria.append(key)

with open('list_of_species.txt', 'w') as tofile:
    for item in gammaproteobacteria[0:5] + proteobacteria[0:5] + bacteria[0:5]:
        tofile.write(item + '\n')

# while read p ; do scp -r fas31@abacus:/data/genomes/bacteria2016-04-14/$p.embl .; done < list_of_species.txt
# for filename in ../data/embl/outgroups/*.embl; do embl2peptides.py $filename ../data/fasta-protein/outgroups/"$(basename $filename | cut -d. -f1).fasta"; done