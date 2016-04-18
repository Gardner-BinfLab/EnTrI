iipath = '../results/check-biases/without-ends.txt'
operonpath = '../operons/SL1344-cropped.opr'
locustag = 'SL1344'
iidict = {}
with open(iipath, 'r') as iifile:
    for iiline in iifile:
        if iiline.startswith(locustag):
            iicells = iiline.split()
            iidict[iicells[0]] = float(iicells[1])
print len(iidict)
print sum(0 <= iidict[x] < 0.2 for x in iidict.keys())


# operonid = ''
# operongenes = []
# operonstrand = ''
operongenesnum = 0
operonessentialsnum = 0
with open(operonpath, 'r') as operonfile:
    for operonline in operonfile:
        operoncells = operonline.split()
        if operoncells[2].startswith(locustag):
            operongenesnum += 1
            if operoncells[2] in iidict.keys() and 0 <= iidict[operoncells[2]] < 0.2:
                operonessentialsnum += 1
print operongenesnum
print operonessentialsnum
