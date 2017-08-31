from os import listdir, remove, system, makedirs
from shutil import copyfile
deg = '/home/fatemeh/EnTrI/results/deg/define-core-accessory-hieranoid-fitch/core-essential-genomes/hhmake/'
endo = '/home/fatemeh/EnTrI/results/endosymbionts/define-core-accessory-hieranoid/core-essential-genomes/hhmake/'
entero = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials/hhmake/'
outdir = '/home/fatemeh/EnTrI/results/venn-entero-deg-endo/'

list_of_degs = listdir(deg)
list_of_endos = listdir(endo)
list_of_enteros = listdir(entero)

# list_of_degs.remove('degfam.hmm')
# list_of_endos.remove('endofam.hmm')
# list_of_enteros.remove('enterofam.hmm')

temppath = deg + 'temp.txt'
for filename in list_of_degs:
    degpath = deg + filename
    copyfile(degpath, temppath)
    with open(temppath, 'r') as fromfile:
        with open(degpath, 'w') as tofile:
            for line in fromfile:
                cells = line.split()
                if len(cells) and cells[0] == 'NAME':
                    cells[1] = filename.split('.')[0]
                    line = '  '.join(cells) + '\n'
                tofile.write(line)
remove(temppath)
system('cat ' + deg + '* > ' + deg + 'degfam.hmm')

temppath = endo + 'temp.txt'
for filename in list_of_endos:
    endopath = endo + filename
    copyfile(endopath, temppath)
    with open(temppath, 'r') as fromfile:
        with open(endopath, 'w') as tofile:
            for line in fromfile:
                cells = line.split()
                if len(cells) and cells[0] == 'NAME':
                    cells[1] = filename.split('.')[0]
                    line = '  '.join(cells) + '\n'
                tofile.write(line)
remove(temppath)
system('cat ' + endo + '* > ' + endo + 'endofam.hmm')

temppath = entero + 'temp.txt'
for filename in list_of_enteros:
    enteropath = entero + filename
    copyfile(enteropath, temppath)
    with open(temppath, 'r') as fromfile:
        with open(enteropath, 'w') as tofile:
            for line in fromfile:
                cells = line.split()
                if len(cells) and cells[0] == 'NAME':
                    cells[1] = filename.split('.')[0]
                    line = '  '.join(cells) + '\n'
                tofile.write(line)
remove(temppath)
system('cat ' + entero + '* > ' + entero + 'enterofam.hmm')

enterodeg = outdir + 'enterodeg/'
makedirs(enterodeg)
for filename in list_of_enteros:
    system('hhsearch -i ' + entero + filename + ' -d ' + deg + 'degfam.hmm -o ' + enterodeg + filename)

enteroendo = outdir + 'enteroendo/'
makedirs(enteroendo)
for filename in list_of_enteros:
    system('hhsearch -i ' + entero + filename + ' -d ' + endo + 'endofam.hmm -o ' + enteroendo + filename)

endodeg = outdir + 'endodeg/'
makedirs(endodeg)
for filename in list_of_endos:
    system('hhsearch -i ' + endo + filename + ' -d ' + deg + 'degfam.hmm -o ' + endodeg + filename)

enterodeg_dict = {}
enteroendo_dict = {}
endodeg_dict = {}
enterodeg_write = outdir + 'entero1deg1endo0.txt'
enteroendo_write = outdir + 'entero1deg0endo1.txt'
endodeg_write = outdir + 'entero0deg1endo1.txt'
enterodegendo_write = outdir + 'entero1deg1endo1.txt'
entero1deg1endo0_write = outdir + 'entero1deg1endo0.txt'
entero1deg0endo1_write = outdir + 'entero1deg0endo1.txt'
entero0deg1endo1_write = outdir + 'entero0deg1endo1.txt'
entero1deg0endo0_write = outdir + 'entero1deg0endo0.txt'
entero0deg1endo0_write = outdir + 'entero0deg1endo0.txt'
entero0deg0endo1_write = outdir + 'entero0deg0endo1.txt'

for filename in list_of_enteros:
    with open(enterodeg + filename, 'r') as fromfile:
        enteroclust = filename.split('.')[0]
        for line in fromfile:
            cells = line.split()
            if len(cells) > 1 and cells[1].startswith('deg-clust'):
                degclust = cells[1]
                prob = float(cells[2])
                if prob >= 99:
                    if degclust not in enterodeg_dict.keys():
                        enterodeg_dict[degclust] = (enteroclust, prob)
                    elif prob > enterodeg_dict[degclust][1]:
                            enterodeg_dict[degclust] = (enteroclust, prob)
                break

for filename in list_of_enteros:
    with open(enteroendo + filename, 'r') as fromfile:
        enteroclust = filename.split('.')[0]
        for line in fromfile:
            cells = line.split()
            if len(cells) > 1 and cells[1].startswith('endo-clust'):
                endoclust = cells[1]
                prob = float(cells[2])
                if prob >= 99:
                    if endoclust not in enteroendo_dict.keys():
                        enteroendo_dict[endoclust] = (enteroclust, prob)
                    elif prob > enteroendo_dict[endoclust][1]:
                            enteroendo_dict[endoclust] = (enteroclust, prob)
                break

for filename in list_of_endos:
    with open(endodeg + filename, 'r') as fromfile:
        endoclust = filename.split('.')[0]
        for line in fromfile:
            cells = line.split()
            if len(cells) > 1 and cells[1].startswith('deg-clust'):
                degclust = cells[1]
                prob = float(cells[2])
                if prob >= 99:
                    if degclust not in endodeg_dict.keys():
                        endodeg_dict[degclust] = (endoclust, prob)
                    elif prob > endodeg_dict[degclust][1]:
                            endodeg_dict[degclust] = (endoclust, prob)
                break

union_entero = []
union_endo = []
union_deg = []
with open(enterodegendo_write, 'w') as tofile:
    for key in enterodeg_dict.keys():
        if key in endodeg_dict.keys() and endodeg_dict[key][0] in enteroendo_dict.keys() and \
                        enteroendo_dict[endodeg_dict[key][0]][0] == enterodeg_dict[key][0]:
            tofile.write(key + '\t' + endodeg_dict[key][0] + '\t' + enterodeg_dict[key][0] + '\n')
            union_entero.append(enterodeg_dict[key][0])
            union_endo.append(endodeg_dict[key][0])
            union_deg.append(key)

with open(enterodeg_write, 'w') as tofile:
    for key in enterodeg_dict.keys():
        if key not in union_deg:
            tofile.write(key + '\t' + enterodeg_dict[key][0] + '\t' + str(enterodeg_dict[key][1]) + '\n')

with open(enteroendo_write, 'w') as tofile:
    for key in enteroendo_dict.keys():
        if key not in union_endo:
            tofile.write(key + '\t' + enteroendo_dict[key][0] + '\t' + str(enteroendo_dict[key][1]) + '\n')

with open(endodeg_write, 'w') as tofile:
    for key in endodeg_dict.keys():
        if key not in union_deg:
            tofile.write(key + '\t' + endodeg_dict[key][0] + '\t' + str(endodeg_dict[key][1]) + '\n')

list_of_enteros = [item.split('.')[0] for item in list_of_enteros]
# print('enteros: ' + str(len(list_of_enteros)))
list_of_endos = [item.split('.')[0] for item in list_of_endos]
# print('endos: ' + str(len(list_of_endos)))
list_of_degs = [item.split('.')[0] for item in list_of_degs]
# print('degs: ' + str(len(list_of_degs)))

with open(entero1deg0endo0_write, 'w') as tofile:
    for key in list_of_enteros:
        if key not in [item[0] for item in list(enteroendo_dict.values())] and key not in [item[0] for item in list(enterodeg_dict.values())]:
            tofile.write(key + '\n')

with open(entero0deg1endo0_write, 'w') as tofile:
    for key in list_of_degs:
        if key not in endodeg_dict.keys() and key not in enterodeg_dict.keys():
            tofile.write(key + '\n')

with open(entero0deg0endo1_write, 'w') as tofile:
    for key in list_of_endos:
        if key not in [item[0] for item in list(endodeg_dict.values())] and key not in enteroendo_dict.keys():
            tofile.write(key + '\n')