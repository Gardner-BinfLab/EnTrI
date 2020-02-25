import pandas as pd
import statistics
from statistics import mode
from os import listdir, path, mkdir
from shutil import rmtree

def makedir(dirname):
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

enteroh = '../results/hieranoid/clusters.txt'
degh = '../results/deg/clusters.txt'
endoh = '../results/endosymbionts/clusters.txt'

enterog = '../results/eggnog-mapper/seqdb.fasta.emapper.annotations'
degg = '../results/deg/fasta-proteins/seqdb.fasta.emapper.annotations'
endog = '../results/endosymbionts/fasta-proteins/seqdb.fasta.emapper.annotations'

outdir = '../results/giant-tab/'
makedir(outdir)
enteroa = outdir + 'entero_annotated_clusters.txt'
dega = outdir + 'deg_annotated_clusters.txt'
endoa = outdir + 'endo_annotated_clusters.txt'

eggnog = pd.read_csv(enterog, sep='\t', header=None).iloc[:,[0, 4, 8]]
eggnog.columns = ['locus','gene', 'bactNOG']
eggnog.set_index('locus',inplace=True)
eggnog.loc[:,'bactNOG'] = eggnog.loc[:,'bactNOG'].str.extract(r'.*(\w{5}@bactNOG).*', expand=False)
eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].isnull(),'gene'] = eggnog.loc[eggnog.loc[:,'gene'].isnull(),'bactNOG']

with open(enteroh, 'r') as fromfile:
    with open(enteroa, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            cells = [i for i in cells if not i.startswith('ENC')]
            if len(cells):
                genes = eggnog['gene'][eggnog.index.isin(cells)].values
                bactnogs = eggnog['bactNOG'][eggnog.index.isin(cells)].values
                if len(genes) or len(bactnogs):
                    try:
                        gene_name = mode(genes)
                    except:
                        gene_name = min(genes, key=len) # to prefer gene names to bactNOG IDs
                    try:
                        nog_name = mode(bactnogs)
                    except:
                        nog_name = bactnogs[0] # randomly select one
                    if nog_name or gene_name:
                        tofile.write(str(nog_name) + '\t' + str(gene_name))
                        for cell in cells:
                            tofile.write('\t' + cell)
                        tofile.write('\n')

eggnog = pd.read_csv(degg, sep='\t', header=None).iloc[:,[0, 4, 9]]
eggnog.columns = ['locus','gene', 'bactNOG']
eggnog.set_index('locus',inplace=True)
eggnog.loc[:,'bactNOG'] = eggnog.loc[:,'bactNOG'].str.extract(r'.*(\w{5}@bactNOG).*', expand=False)
eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].isnull(),'gene'] = eggnog.loc[eggnog.loc[:,'gene'].isnull(),'bactNOG']

with open(degh, 'r') as fromfile:
    with open(dega, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            cells = [i for i in cells if not i.startswith('DEG1019')]
            cells = [i for i in cells if not i.startswith('DEG1032')]
            cells = [i for i in cells if not i.startswith('DEG1033')]
            cells = [i for i in cells if not i.startswith('DEG1048')]
            if len(cells):
                genes = eggnog['gene'][eggnog.index.isin(cells)].values
                bactnogs = eggnog['bactNOG'][eggnog.index.isin(cells)].values
                if len(genes) or len(bactnogs):
                    try:
                        gene_name = mode(genes)
                    except:
                        gene_name = min(genes, key=len) # to prefer gene names to bactNOG IDs
                    try:
                        nog_name = mode(bactnogs)
                    except:
                        nog_name = bactnogs[0] # randomly select one
                    if nog_name or gene_name:
                        tofile.write(str(nog_name) + '\t' + str(gene_name))
                        for cell in cells:
                            tofile.write('\t' + cell)
                        tofile.write('\n')

eggnog = pd.read_csv(endog, sep='\t', header=None).iloc[:,[0, 4, 9]]
eggnog.columns = ['locus','gene', 'bactNOG']
eggnog.set_index('locus',inplace=True)
eggnog.loc[:,'bactNOG'] = eggnog.loc[:,'bactNOG'].str.extract(r'.*(\w{5}@bactNOG).*', expand=False)
eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()>4,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'gene'] = eggnog.loc[eggnog.loc[:,'gene'].str.len()<3,'bactNOG']
eggnog.loc[eggnog.loc[:,'gene'].isnull(),'gene'] = eggnog.loc[eggnog.loc[:,'gene'].isnull(),'bactNOG']

with open(endoh, 'r') as fromfile:
    with open(endoa, 'w') as tofile:
        for line in fromfile:
            cells = line.split()
            if len(cells):
                genes = eggnog['gene'][eggnog.index.isin(cells)].values
                bactnogs = eggnog['bactNOG'][eggnog.index.isin(cells)].values
                if len(genes) or len(bactnogs):
                    try:
                        gene_name = mode(genes)
                    except:
                        try:
                            gene_name = min(genes, key=len) # to prefer gene names to bactNOG IDs
                        except:
                            gene_name = genes[0]
                    try:
                        nog_name = mode(bactnogs)
                    except:
                        nog_name = bactnogs[0] # randomly select one
                    if nog_name or gene_name:
                        tofile.write(str(nog_name) + '\t' + str(gene_name))
                        for cell in cells:
                            tofile.write('\t' + cell)
                        tofile.write('\n')