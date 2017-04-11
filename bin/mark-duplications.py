from os import listdir
from re import match

interesting_genes = '../results/interesting_genes/sometimes-essential.tsv'
no_duplications = '../results/interesting_genes/sometimes-essential-marked-dup.tsv'
names = ['ENC', 'BN373', 'ERS227112', 'ROD', 'SL1344', 'SL3261', 'STMMW', 'STM', 'SEN', 't', 'EC958', 'NCTC13441',
         'CS17', 'ETEC', 'BW25113', 'b']
clusters = '../results/homclust/EFam-clusters/'
list_of_files = listdir(clusters)
with open(interesting_genes, 'r') as fromfile:
    with open(no_duplications, 'w') as tofile:
        lines = fromfile.readlines()
        tofile.write(lines[0].rstrip('\n') + '\tENC_no\tBN373_no\tERS227112_no\tROD_no\tSL1344_no\tSL3261_no\tSTMMW_no\tSTM_no\t'
                                + 'SEN_no\tt_no\tEC958_no\tNCTC13441_no\tCS17_no\tETEC_no\tBW25113_no\tb_no\n')
        for line in lines[1:len(lines)]:
            cells = line.split('\t')
            genes = cells[18:34]
            essentialities = {'ENC': 0, 'BN373': 0, 'ERS227112': 0, 'ROD': 0, 'SL1344': 0, 'SL3261': 0, 'STMMW': 0,
                              'STM': 0, 'SEN': 0, 't': 0, 'EC958': 0, 'NCTC13441': 0, 'CS17': 0, 'ETEC': 0,
                              'BW25113': 0, 'b': 0}
            for i in range(2, 18):
                if cells[i] != 'X' and float(cells[i]) >= 1.644854:
                    essentialities[names[i - 2]] += 1
            i = 0
            first_gene = genes[i]
            while first_gene == 'X':
                i += 1
                first_gene = genes[i]
            list_of_clusters = list()
            for filename in list_of_files:
                with open(clusters + filename, 'r') as clusterfile:
                    for row in clusterfile:
                        if first_gene in row:
                            list_of_clusters.append(filename)
                            break
            genomes = {'ENC': 0, 'BN373': 0, 'ERS227112': 0, 'ROD': 0, 'SL1344': 0, 'SL3261': 0, 'STMMW': 0, 'STM': 0,
                       'SEN': 0, 't': 0, 'EC958': 0, 'NCTC13441': 0, 'CS17': 0, 'ETEC': 0, 'BW25113': 0, 'b': 0}
            for filename in list_of_clusters:
                with open(clusters + filename, 'r') as clusterfile:
                    for row in clusterfile:
                        items = row.split()
                        genome = match('([a-zA-Z0-9]+_+|[a-zA-Z]+)[a-zA-Z0-9]+', items[1]).group(1).strip('_')
                        genomes[genome] += 1
            # i = 0
            # item = names[i]
            # while genomes[item] > 1 or (essentialities[item] and genomes[item]):
            #     i += 1
            #     if i < len(names):
            #         item = names[i]
            #     else:
            #         break
            # # if i != len(names):
            # #     tofile.write(line)
            # if i == len(names):
            #     line = str.replace(line, cells[0], cells[0]+'_dup', 1)
            line = line.rstrip('\n') + '\t' + str(genomes['ENC']) + '\t' + str(genomes['BN373']) + '\t' + str(genomes['ERS227112'])\
                   + '\t' + str(genomes['ROD']) + '\t' + str(genomes['SL1344']) + '\t' + str(genomes['SL3261']) + '\t'\
                   + str(genomes['STMMW']) + '\t' + str(genomes['STM']) + '\t' + str(genomes['SEN']) + '\t' +\
                   str(genomes['t']) + '\t' + str(genomes['EC958']) + '\t' + str(genomes['NCTC13441']) + '\t' +\
                   str(genomes['CS17']) + '\t' + str(genomes['ETEC']) + '\t' + str(genomes['BW25113']) + '\t' +\
                   str(genomes['b']) + '\n'
            tofile.write(line)