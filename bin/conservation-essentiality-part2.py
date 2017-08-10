from os import listdir
# def read_gene_essentiality(indir):
#     list_of_files = listdir(indir)
#     iidict = {}
#     for filename in list_of_files:
#         with open(indir+'/'+filename) as from_file:
#             for line in from_file:
#                 cells = line.split()
#                 iidict[cells[0]] = float(cells[1])
#     return iidict

# iipath = '/home/fatemeh/EnTrI/results/biases/dbscan'
# insertion_indices = read_gene_essentiality(iipath)
esspath = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials' #
esslist = listdir(esspath) #
esslist.remove('profiles') #
esslist.remove('alignments') #
clusters = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-cores/'
writepath = '/home/fatemeh/EnTrI/results/homologs.essentiality'
# num_species = 14

list_of_files = listdir(clusters)
with open(writepath, 'w') as tofile:
    for filename in list_of_files:
        if filename.startswith('clust'):
            if filename in esslist: #
                essentiality = 1 #
            else: #
                essentiality = 0 #
            # with open(clusters + filename, 'r') as fromfile:
            #     essentiality = 0
            #     count = 0
            #     for line in fromfile:
            #         if line.startswith('>'):
            #             gene = line.strip()[1:]
            #             if gene in insertion_indices.keys():
            #                 essentiality = 1
            #                 essentiality += insertion_indices[gene]
            #                 count += 1
            # essentiality /= count
            tofile.write(filename + '\t' + str(essentiality) + '\n')
