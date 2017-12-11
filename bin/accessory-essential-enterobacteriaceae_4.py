genes = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap-k12.tsv'
kegg = '/home/fatemeh/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat'
outpath = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap-kegg.tsv'

kegg_dict = {}
with open(kegg, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[2] = cells[2][1:-2]
        if cells[0] != 'gene_id' and cells[0] not in kegg_dict.keys():
            kegg_dict[cells[0]] = cells[2]
        elif cells[0] != 'gene_id':
            kegg_dict[cells[0]] += ' / ' + cells[2]

with open(genes, 'r') as fromfile:
    with open(outpath, 'w') as tofile:
        for line in fromfile:
            gene = line.rstrip()
            if gene in kegg_dict.keys():
                tofile.write(gene + '\t' + kegg_dict[gene] + '\n')
            else:
                tofile.write(gene + '\t' + '\n')