from scipy.stats import fisher_exact

hmmscanpath = '../results/beneficialloss-plasmid/hmmscan.domtblout'
hmmscan = []
with open(hmmscanpath, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('#'):
            cells = line.split()
            hmmscan.append(cells[0])
hmmscan = list(set(hmmscan))

benlosspath = '../results/beneficialloss-plasmid/beneficial_losses.txt'
benloss = []
with open(benlosspath, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('.')
        benloss.append(cells[0])

otherspath = '../results/beneficialloss-plasmid/others.txt'
others = []
with open(otherspath, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('.')
        others.append(cells[0])

benloss_copied = 0
benloss_not = 0
others_copied = 0
others_not = 0

for item in benloss:
    if item in hmmscan:
        benloss_copied += 1
    else:
        benloss_not += 1

for item in others:
    if item in hmmscan:
        others_copied += 1
    else:
        others_not += 1

oddsratio, pvalue = fisher_exact([[benloss_copied, benloss_not], [others_copied, others_not]])
print(pvalue)
print(benloss_copied)
print(benloss_not)
print(others_copied)
print(others_not)