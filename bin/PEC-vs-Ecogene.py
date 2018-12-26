ecopath = '../results/ecogene-k12.txt'
pecpath = '../results/PEC.txt'
eco = set()
with open(ecopath, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        eco.add(cells[0])

pec = set()
with open(pecpath, 'r') as fromfile:
    fromfile.readline()
    for line in fromfile:
        cells = line.split()
        newcells = cells[3].split(',')
        i = 0
        while not(newcells[i].startswith('b')):
            i += 1
        pec.add(newcells[i])

pde = pec - eco
dpe = eco - pec
print(len(pde))
print(len(dpe))
intersect = eco.intersection(pec)
print(len(intersect))