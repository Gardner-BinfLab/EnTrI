eggnog = '../results/define-core-accessory-hieranoid-fitch/genus-seqdb.fasta.emapper.annotations'
gene_dict = {}
with open(eggnog, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        id = cells[0]
        gene = cells[4]
        if gene == '':
            all = cells[9].split(',')
            for item in all:
                pre = item.split('@')[0]
                post = item.split('@')[1]
                if post == 'NOG':
                    gene = 'eggNOG_' + pre
                    break
        else:
            first = gene[0:3]
            last = gene[3:]
            gene = first.lower() + last
        if gene not in gene_dict.keys():
            gene_dict[gene] = [0, 0, 0, 0, 0] # 0: Citrobacter, 1:Salmonella, 2: Escherichia, 3: Klebsiella, 4: Enterobacter
        if id.startswith('ROD'):
            gene_dict[gene][0] = 1
        elif id.startswith('t'):
            gene_dict[gene][1] = 1
        elif id.startswith('b'):
            gene_dict[gene][2] = 1
        elif id.startswith('BN373'):
            gene_dict[gene][3] = 1
        elif id.startswith('ENC'):
            gene_dict[gene][4] = 1

for item in gene_dict.keys():
    gene_dict[item] = tuple(gene_dict[item])

sorted_dict = sorted(gene_dict, key=gene_dict.get, reverse=True)
writepath = '../results/define-core-accessory-hieranoid-fitch/heatmap.tsv'
with open(writepath, 'w') as tofile:
    tofile.write('gene\tCitrobacter\tSalmonella\tEscherichia\tKlebsiella\tEnterobacter\n')
    for item in sorted_dict:
        tofile.write(item)
        for i in gene_dict[item]:
            tofile.write('\t' + str(i))
        tofile.write('\n')