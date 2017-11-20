clusters = '../results/venn-entero-deg-endo/final-clusters-plus-single-groups.txt'
eggnog = '../results/venn-entero-deg-endo/final-clusters-plus-single-groups.fasta.emapper.annotations'
heatmap = '../results/venn-entero-deg-endo/heatmap.tsv'
gene_dict = {}

with open(heatmap, 'w') as tofile:
    with open(clusters, 'r') as clustfile:
        with open(eggnog, 'r') as eggnogfile:
            tofile.write('gene\tEnterobacteriaceae\tEndosymbiont\tGammaproteobacteria\tProteobacteria\tBacteria\n')
            for line in clustfile:
                clades = line.split()
                clades = clades[1:]
                for i in range(0, len(clades)):
                    clades[i] = clades[i].split('-')[0]
                eggnogline = eggnogfile.readline()
                cells = eggnogline.split('\t')
                first = cells[4][0:3]
                last = cells[4][3:]
                gene = first.lower() + last
                if gene == '':
                    gene = cells[9]. split('@')[0]
                # tofile.write(gene + '\t' + str(int('entero' in clades)) + '\t' + str(int('endo' in clades)) + '\t' +
                #              str(int('gamma' in clades)) + '\t' + str(int('prot' in clades)) + '\t' +
                #              str(int('deg' in clades)) + '\n')
                gene_dict[gene] = (int('entero' in clades), int('endo' in clades), int('gamma' in clades),
                                   int('prot' in clades), int('deg' in clades))
            for gene in sorted(gene_dict, key=gene_dict.get, reverse=True):
                tofile.write(gene)
                for item in gene_dict[gene]:
                    tofile.write('\t' + str(item))
                tofile.write('\n')