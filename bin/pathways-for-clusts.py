from os import listdir
# from Bio import SeqIO

def read_fasta_sequences(filepath):
    fasta_dict = {}
    with open(filepath, 'r') as fasta_file:
        for line in fasta_file:
            cells = line.split('\t')
            fasta_dict[cells[0]] = cells[4]+'\t'+cells[-1]
    return fasta_dict

genes = {}
coreclusts = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-cores/'
list_of_files = listdir(coreclusts)
list_of_files = [item for item in list_of_files if item.startswith('clust')]
for filename in list_of_files:
    with open(coreclusts + filename, 'r') as fromfile:
        gene = ''
        for line in fromfile:
            if line.startswith('>b'):
                gene = line.strip()[1:]
                break
        if gene != '':
            genes[gene] = filename

writepath = '/home/fatemeh/EnTrI/results/homologs.path'
keggpath = '/home/fatemeh/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat'
with open(keggpath, 'r') as fromfile:
    with open(writepath, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            gene = cells[0].strip('"')
            if gene == 'gene_id':
                tofile.write(line)
            elif gene in genes.keys():
                tofile.write('"'+genes[gene]+'"\t'+cells[1]+'\t'+' '.join(cells[2:])+'\n')

# writepath = '/home/fatemeh/endosymbionts/essential-not-in-endosymbionts.txt'
# with open(writepath, 'w') as tofile:
#     for item in newclusts:
#         gene = ''
#         with open(clusters + item, 'r') as fromfile:
#             for line in fromfile:
#                 if line.startswith('>b'):
#                     gene = line.strip()[1:]
#                     break
#         tofile.write(gene+'\t'+fastas[gene])