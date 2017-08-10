from os import listdir
# from Bio import SeqIO

def read_eggnog(filepath):
    fasta_dict = {}
    with open(filepath, 'r') as fasta_file:
        for line in fasta_file:
            cells = line.split('\t')
            if cells[4] == '':
                cells[4] = '-'
            fasta_dict[cells[0]] = cells[4]+'\t'+cells[-1][0:-1]
    return fasta_dict


def read_kegg(filepath):
    path_dict = {}
    with open(filepath, 'r') as fromfile:
        for line in fromfile:
            cells = line.split('\t')
            for i in range(0, len(cells)):
                cells[i] = cells[i].strip('"')
            if cells[0] in path_dict.keys():
                addition = ' / ' + cells[2][0:-2]
                path_dict[cells[0]] += addition
            else:
                path_dict[cells[0]] = cells[2][0:-2]
    return path_dict

clusters = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials/'
list_of_clusts = listdir(clusters)
list_of_clusts.remove('alignments')
list_of_clusts.remove('profiles')

writepath = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials/clusters.info'
coreessfasta = '/home/fatemeh/EnTrI/results/eggnog-mapper/U00096.fasta.emapper.annotations'
eggnog = read_eggnog(coreessfasta)
keggpath = '/home/fatemeh/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat'
kegg = read_kegg(keggpath)
print(kegg)
with open(writepath, 'w') as tofile:
    for key in list_of_clusts:
        gene = ''
        with open(clusters + key, 'r') as fromfile:
            for line in fromfile:
                if line.startswith('>b'):
                    gene = line.strip()[1:]
                    break
        try:
            tofile.write(key+'\t'+eggnog[gene] + '\t' + kegg[gene] + '\n')
        except:
            try:
                tofile.write(key + '\t' + eggnog[gene] + '\t' + '-' + '\n')
            except:
                try:
                    tofile.write(key + '\t' + '-\t-' + '\t' + kegg[gene] + '\n')
                except:
                    tofile.write(key + '\t' + '-\t-' + '\t' + '-' + '\n')