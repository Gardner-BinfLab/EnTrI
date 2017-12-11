heatmap = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap.tsv'
eggnog = '/home/fatemeh/EnTrI/results/eggnog-mapper/U00096.fasta.emapper.annotations'
outpath = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap-k12.tsv'

gene_dict = {}
with open(eggnog, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        if cells[4] != '':
            gene_dict[cells[4]] = cells[0]

with open(heatmap, 'r') as fromfile:
    with open(outpath, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            if cells[0] != 'gene':
                genename = cells[0].upper()
                if genename in gene_dict.keys():
                    tofile.write(gene_dict[genename] + '\n')
                else:
                    tofile.write('NA' + '\n')
