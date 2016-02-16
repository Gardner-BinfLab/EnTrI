#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

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

seqdb = '/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/merge-clust-plot'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory'
esscoredir = outdir + '/core-essential-genomes'
sesscoredir = outdir + '/core-sometimes-essential-genomes'
nesscoredir = outdir + '/core-never-essential-genomes'
essaccessorydir = outdir + '/accessory-essential-genomes'
sessaccessorydir = outdir + '/accessory-sometimes-essential-genomes'
nessaccessorydir = outdir + '/accessory-never-essential-genomes'

salesscoredir = outdir + '/salmonella-core-essential-genomes'
salsesscoredir = outdir + '/salmonella-core-sometimes-essential-genomes'
salnesscoredir = outdir + '/salmonella-core-never-essential-genomes'
salessaccessorydir = outdir + '/salmonella-accessory-essential-genomes'
salsessaccessorydir = outdir + '/salmonella-accessory-sometimes-essential-genomes'
salnessaccessorydir = outdir + '/salmonella-accessory-never-essential-genomes'

ecoesscoredir = outdir + '/ecoli-core-essential-genomes'
ecosesscoredir = outdir + '/ecoli-core-sometimes-essential-genomes'
econesscoredir = outdir + '/ecoli-core-never-essential-genomes'
ecoessaccessorydir = outdir + '/ecoli-accessory-essential-genomes'
ecosessaccessorydir = outdir + '/ecoli-accessory-sometimes-essential-genomes'
econessaccessorydir = outdir + '/ecoli-accessory-never-essential-genomes'

kleesscoredir = outdir + '/klebsiella-core-essential-genomes'
klesesscoredir = outdir + '/klebsiella-core-sometimes-essential-genomes'
klenesscoredir = outdir + '/klebsiella-core-never-essential-genomes'
kleessaccessorydir = outdir + '/klebsiella-accessory-essential-genomes'
klesessaccessorydir = outdir + '/klebsiella-accessory-sometimes-essential-genomes'
klenessaccessorydir = outdir + '/klebsiella-accessory-never-essential-genomes'

citesscoredir = outdir + '/citrobacter-core-essential-genomes'
citsesscoredir = outdir + '/citrobacter-core-sometimes-essential-genomes'
citnesscoredir = outdir + '/citrobacter-core-never-essential-genomes'
citessaccessorydir = outdir + '/citrobacter-accessory-essential-genomes'
citsessaccessorydir = outdir + '/citrobacter-accessory-sometimes-essential-genomes'
citnessaccessorydir = outdir + '/citrobacter-accessory-never-essential-genomes'

entesscoredir = outdir + '/enterobacter-core-essential-genomes'
entsesscoredir = outdir + '/enterobacter-core-sometimes-essential-genomes'
entnesscoredir = outdir + '/enterobacter-core-never-essential-genomes'
entessaccessorydir = outdir + '/enterobacter-accessory-essential-genomes'
entsessaccessorydir = outdir + '/enterobacter-accessory-sometimes-essential-genomes'
entnessaccessorydir = outdir + '/enterobacter-accessory-never-essential-genomes'

seesscoredir = outdir + '/salmonella-ecoli-core-essential-genomes'
sesesscoredir = outdir + '/salmonella-ecoli-core-sometimes-essential-genomes'
senesscoredir = outdir + '/salmonella-ecoli-core-never-essential-genomes'
seessaccessorydir = outdir + '/salmonella-ecoli-accessory-essential-genomes'
sesessaccessorydir = outdir + '/salmonella-ecoli-accessory-sometimes-essential-genomes'
senessaccessorydir = outdir + '/salmonella-ecoli-accessory-never-essential-genomes'

secesscoredir = outdir + '/salmonella-ecoli-citrobacter-core-essential-genomes'
secsesscoredir = outdir + '/salmonella-ecoli-citrobacter-core-sometimes-essential-genomes'
secnesscoredir = outdir + '/salmonella-ecoli-citrobacter-core-never-essential-genomes'
secessaccessorydir = outdir + '/salmonella-ecoli-citrobacter-accessory-essential-genomes'
secsessaccessorydir = outdir + '/salmonella-ecoli-citrobacter-accessory-sometimes-essential-genomes'
secnessaccessorydir = outdir + '/salmonella-ecoli-citrobacter-accessory-never-essential-genomes'

keesscoredir = outdir + '/klebsiella-enterobacter-core-essential-genomes'
kesesscoredir = outdir + '/klebsiella-enterobacter-core-sometimes-essential-genomes'
kenesscoredir = outdir + '/klebsiella-enterobacter-core-never-essential-genomes'
keessaccessorydir = outdir + '/klebsiella-enterobacter-accessory-essential-genomes'
kesessaccessorydir = outdir + '/klebsiella-enterobacter-accessory-sometimes-essential-genomes'
kenessaccessorydir = outdir + '/klebsiella-enterobacter-accessory-never-essential-genomes'

makedir(outdir)
makedir(esscoredir)
makedir(sesscoredir)
makedir(nesscoredir)
makedir(essaccessorydir)
makedir(sessaccessorydir)
makedir(nessaccessorydir)

makedir(salesscoredir)
makedir(salsesscoredir)
makedir(salnesscoredir)
makedir(salessaccessorydir)
makedir(salsessaccessorydir)
makedir(salnessaccessorydir)

makedir(ecoesscoredir)
makedir(ecosesscoredir)
makedir(econesscoredir)
makedir(ecoessaccessorydir)
makedir(ecosessaccessorydir)
makedir(econessaccessorydir)

makedir(kleesscoredir)
makedir(klesesscoredir)
makedir(klenesscoredir)
makedir(kleessaccessorydir)
makedir(klesessaccessorydir)
makedir(klenessaccessorydir)

makedir(citesscoredir)
makedir(citsesscoredir)
makedir(citnesscoredir)
makedir(citessaccessorydir)
makedir(citsessaccessorydir)
makedir(citnessaccessorydir)

makedir(entesscoredir)
makedir(entsesscoredir)
makedir(entnesscoredir)
makedir(entessaccessorydir)
makedir(entsessaccessorydir)
makedir(entnessaccessorydir)

makedir(seesscoredir)
makedir(sesesscoredir)
makedir(senesscoredir)
makedir(seessaccessorydir)
makedir(sesessaccessorydir)
makedir(senessaccessorydir)

makedir(secesscoredir)
makedir(secsesscoredir)
makedir(secnesscoredir)
makedir(secessaccessorydir)
makedir(secsessaccessorydir)
makedir(secnessaccessorydir)

makedir(keesscoredir)
makedir(kesesscoredir)
makedir(kenesscoredir)
makedir(keessaccessorydir)
makedir(kesessaccessorydir)
makedir(kenessaccessorydir)

esscoregenes = []
sesscoregenes = []
nesscoregenes = []
essaccessorygenes = []
sessaccessorygenes = []
nessaccessorygenes = []

salesscoregenes = []
salsesscoregenes = []
salnesscoregenes = []
salessaccessorygenes = []
salsessaccessorygenes = []
salnessaccessorygenes = []

ecoesscoregenes = []
ecosesscoregenes = []
econesscoregenes = []
ecoessaccessorygenes = []
ecosessaccessorygenes = []
econessaccessorygenes = []

kleesscoregenes = []
klesesscoregenes = []
klenesscoregenes = []
kleessaccessorygenes = []
klesessaccessorygenes = []
klenessaccessorygenes = []

citesscoregenes = []
citsesscoregenes = []
citnesscoregenes = []
citessaccessorygenes = []
citsessaccessorygenes = []
citnessaccessorygenes = []

entesscoregenes = []
entsesscoregenes = []
entnesscoregenes = []
entessaccessorygenes = []
entsessaccessorygenes = []
entnessaccessorygenes = []

seesscoregenes = []
sesesscoregenes = []
senesscoregenes = []
seessaccessorygenes = []
sesessaccessorygenes = []
senessaccessorygenes = []

secesscoregenes = []
secsesscoregenes = []
secnesscoregenes = []
secessaccessorygenes = []
secsessaccessorygenes = []
secnessaccessorygenes = []

keesscoregenes = []
kesesscoregenes = []
kenesscoregenes = []
keessaccessorygenes = []
kesessaccessorygenes = []
kenessaccessorygenes = []

species_names = ["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t", "b"]
num_species = len(species_names)
gene_dict = {species_names[i]: 0 for i in range(num_species)}
essentiality_dict = {species_names[i]: 0 for i in range(num_species-1)}

sal_species_names = ["SEN", "SL1344", "STM", "STMMW", "t"]
sal_num_species = len(sal_species_names)
sal_gene_dict = {sal_species_names[i]: 0 for i in range(sal_num_species)}
sal_essentiality_dict = {sal_species_names[i]: 0 for i in range(sal_num_species)}

eco_species_names = ["CS17", "ETEC", "NCTC13441"]
eco_num_species = len(eco_species_names)
eco_gene_dict = {eco_species_names[i]: 0 for i in range(eco_num_species)}
eco_essentiality_dict = {eco_species_names[i]: 0 for i in range(eco_num_species)}

kle_species_names = ["ERS227112", "BN373"]
kle_num_species = len(kle_species_names)
kle_gene_dict = {kle_species_names[i]: 0 for i in range(kle_num_species)}
kle_essentiality_dict = {kle_species_names[i]: 0 for i in range(kle_num_species)}

cit_species_names = ["ROD"]
cit_num_species = len(cit_species_names)
cit_gene_dict = {cit_species_names[i]: 0 for i in range(cit_num_species)}
cit_essentiality_dict = {cit_species_names[i]: 0 for i in range(cit_num_species)}

ent_species_names = ["ENC"]
ent_num_species = len(ent_species_names)
ent_gene_dict = {ent_species_names[i]: 0 for i in range(ent_num_species)}
ent_essentiality_dict = {ent_species_names[i]: 0 for i in range(ent_num_species)}

se_species_names = ["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441"]
se_num_species = len(se_species_names)
se_gene_dict = {se_species_names[i]: 0 for i in range(se_num_species)}
se_essentiality_dict = {se_species_names[i]: 0 for i in range(se_num_species)}

sec_species_names = ["SEN", "SL1344", "STM", "STMMW", "t", "CS17", "ETEC", "NCTC13441", "ROD"]
sec_num_species = len(sec_species_names)
sec_gene_dict = {sec_species_names[i]: 0 for i in range(sec_num_species)}
sec_essentiality_dict = {sec_species_names[i]: 0 for i in range(sec_num_species)}

ke_species_names = ["ERS227112", "BN373", "ENC"]
ke_num_species = len(ke_species_names)
ke_gene_dict = {ke_species_names[i]: 0 for i in range(ke_num_species)}
ke_essentiality_dict = {ke_species_names[i]: 0 for i in range(ke_num_species)}

list_of_files = listdir(clusters)
for filename in list_of_files:
    for key in gene_dict.keys():
        gene_dict[key] = 0
        if key in essentiality_dict.keys():
            essentiality_dict[key] = 0

    for key in sal_gene_dict.keys():
        sal_gene_dict[key] = 0
        sal_essentiality_dict[key] = 0

    for key in eco_gene_dict.keys():
        eco_gene_dict[key] = 0
        eco_essentiality_dict[key] = 0

    for key in kle_gene_dict.keys():
        kle_gene_dict[key] = 0
        kle_essentiality_dict[key] = 0

    for key in cit_gene_dict.keys():
        cit_gene_dict[key] = 0
        cit_essentiality_dict[key] = 0

    for key in ent_gene_dict.keys():
        ent_gene_dict[key] = 0
        ent_essentiality_dict[key] = 0

    for key in se_gene_dict.keys():
        se_gene_dict[key] = 0
        se_essentiality_dict[key] = 0

    for key in sec_gene_dict.keys():
        sec_gene_dict[key] = 0
        sec_essentiality_dict[key] = 0

    for key in ke_gene_dict.keys():
        ke_gene_dict[key] = 0
        ke_essentiality_dict[key] = 0

    list_of_genes = []
    sal_list_of_genes = []
    eco_list_of_genes = []
    kle_list_of_genes = []
    cit_list_of_genes = []
    ent_list_of_genes = []
    se_list_of_genes = []
    sec_list_of_genes = []
    ke_list_of_genes = []

    with open('{0}/{1}'.format(clusters, filename)) as from_file:
        for line in from_file:
            cells = line.split()
            match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
            if not match_result:
                match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
            if match_result:
                name = match_result.group(1)

                if name in gene_dict.keys():
                    list_of_genes.append(cells[1])
                    gene_dict[name] = 1
                    if name in essentiality_dict.keys() and float(cells[4]) < 0.2:
                        essentiality_dict[name] = 1

                if name in sal_gene_dict.keys():
                    sal_list_of_genes.append(cells[1])
                    sal_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        sal_essentiality_dict[name] = 1

                if name in eco_gene_dict.keys():
                    eco_list_of_genes.append(cells[1])
                    eco_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        eco_essentiality_dict[name] = 1

                if name in kle_gene_dict.keys():
                    kle_list_of_genes.append(cells[1])
                    kle_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        kle_essentiality_dict[name] = 1

                if name in cit_gene_dict.keys():
                    cit_list_of_genes.append(cells[1])
                    cit_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        cit_essentiality_dict[name] = 1

                if name in ent_gene_dict.keys():
                    ent_list_of_genes.append(cells[1])
                    ent_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        ent_essentiality_dict[name] = 1

                if name in se_gene_dict.keys():
                    se_list_of_genes.append(cells[1])
                    se_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        se_essentiality_dict[name] = 1

                if name in sec_gene_dict.keys():
                    sec_list_of_genes.append(cells[1])
                    sec_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        sec_essentiality_dict[name] = 1

                if name in ke_gene_dict.keys():
                    ke_list_of_genes.append(cells[1])
                    ke_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        ke_essentiality_dict[name] = 1

    if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
        if essentiality_dict[min(essentiality_dict, key=essentiality_dict.get)] == 1:
            esscoregenes += list_of_genes
        elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
            nesscoregenes += list_of_genes
        else:
            sesscoregenes += list_of_genes
    else:
        if not sum([gene_dict[item] - essentiality_dict[item] for item in essentiality_dict.keys()]):
            essaccessorygenes += list_of_genes
        elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
            nessaccessorygenes += list_of_genes
        else:
            sessaccessorygenes += list_of_genes

    if sal_gene_dict[min(sal_gene_dict, key=sal_gene_dict.get)] == 1:
        if sal_essentiality_dict[min(sal_essentiality_dict, key=sal_essentiality_dict.get)] == 1:
            salesscoregenes += sal_list_of_genes
        elif sal_essentiality_dict[max(sal_essentiality_dict, key=sal_essentiality_dict.get)] == 0:
            salnesscoregenes += sal_list_of_genes
        else:
            salsesscoregenes += sal_list_of_genes
    else:
        if not sum([sal_gene_dict[item] - sal_essentiality_dict[item] for item in sal_essentiality_dict.keys()]):
            salessaccessorygenes += sal_list_of_genes
        elif sal_essentiality_dict[max(sal_essentiality_dict, key=sal_essentiality_dict.get)] == 0:
            salnessaccessorygenes += sal_list_of_genes
        else:
            salsessaccessorygenes += sal_list_of_genes

    if eco_gene_dict[min(eco_gene_dict, key=eco_gene_dict.get)] == 1:
        if eco_essentiality_dict[min(eco_essentiality_dict, key=eco_essentiality_dict.get)] == 1:
            ecoesscoregenes += eco_list_of_genes
        elif eco_essentiality_dict[max(eco_essentiality_dict, key=eco_essentiality_dict.get)] == 0:
            econesscoregenes += eco_list_of_genes
        else:
            ecosesscoregenes += eco_list_of_genes
    else:
        if not sum([eco_gene_dict[item] - eco_essentiality_dict[item] for item in eco_essentiality_dict.keys()]):
            ecoessaccessorygenes += eco_list_of_genes
        elif eco_essentiality_dict[max(eco_essentiality_dict, key=eco_essentiality_dict.get)] == 0:
            econessaccessorygenes += eco_list_of_genes
        else:
            ecosessaccessorygenes += eco_list_of_genes

    if kle_gene_dict[min(kle_gene_dict, key=kle_gene_dict.get)] == 1:
        if kle_essentiality_dict[min(kle_essentiality_dict, key=kle_essentiality_dict.get)] == 1:
            kleesscoregenes += kle_list_of_genes
        elif kle_essentiality_dict[max(kle_essentiality_dict, key=kle_essentiality_dict.get)] == 0:
            klenesscoregenes += kle_list_of_genes
        else:
            klesesscoregenes += kle_list_of_genes
    else:
        if not sum([kle_gene_dict[item] - kle_essentiality_dict[item] for item in kle_essentiality_dict.keys()]):
            kleessaccessorygenes += kle_list_of_genes
        elif kle_essentiality_dict[max(kle_essentiality_dict, key=kle_essentiality_dict.get)] == 0:
            klenessaccessorygenes += kle_list_of_genes
        else:
            klesessaccessorygenes += kle_list_of_genes

    if cit_gene_dict[min(cit_gene_dict, key=cit_gene_dict.get)] == 1:
        if cit_essentiality_dict[min(cit_essentiality_dict, key=cit_essentiality_dict.get)] == 1:
            citesscoregenes += cit_list_of_genes
        elif cit_essentiality_dict[max(cit_essentiality_dict, key=cit_essentiality_dict.get)] == 0:
            citnesscoregenes += cit_list_of_genes
        else:
            citsesscoregenes += cit_list_of_genes
    else:
        if not sum([cit_gene_dict[item] - cit_essentiality_dict[item] for item in cit_essentiality_dict.keys()]):
            citessaccessorygenes += cit_list_of_genes
        elif cit_essentiality_dict[max(cit_essentiality_dict, key=cit_essentiality_dict.get)] == 0:
            citnessaccessorygenes += cit_list_of_genes
        else:
            citsessaccessorygenes += cit_list_of_genes

    if ent_gene_dict[min(ent_gene_dict, key=ent_gene_dict.get)] == 1:
        if ent_essentiality_dict[min(ent_essentiality_dict, key=ent_essentiality_dict.get)] == 1:
            entesscoregenes += ent_list_of_genes
        elif ent_essentiality_dict[max(ent_essentiality_dict, key=ent_essentiality_dict.get)] == 0:
            entnesscoregenes += ent_list_of_genes
        else:
            entsesscoregenes += ent_list_of_genes
    else:
        if not sum([ent_gene_dict[item] - ent_essentiality_dict[item] for item in ent_essentiality_dict.keys()]):
            entessaccessorygenes += ent_list_of_genes
        elif ent_essentiality_dict[max(ent_essentiality_dict, key=ent_essentiality_dict.get)] == 0:
            entnessaccessorygenes += ent_list_of_genes
        else:
            entsessaccessorygenes += ent_list_of_genes

    if se_gene_dict[min(se_gene_dict, key=se_gene_dict.get)] == 1:
        if se_essentiality_dict[min(se_essentiality_dict, key=se_essentiality_dict.get)] == 1:
            seesscoregenes += se_list_of_genes
        elif se_essentiality_dict[max(se_essentiality_dict, key=se_essentiality_dict.get)] == 0:
            senesscoregenes += se_list_of_genes
        else:
            sesesscoregenes += se_list_of_genes
    else:
        if not sum([se_gene_dict[item] - se_essentiality_dict[item] for item in se_essentiality_dict.keys()]):
            seessaccessorygenes += se_list_of_genes
        elif se_essentiality_dict[max(se_essentiality_dict, key=se_essentiality_dict.get)] == 0:
            senessaccessorygenes += se_list_of_genes
        else:
            sesessaccessorygenes += se_list_of_genes

    if sec_gene_dict[min(sec_gene_dict, key=sec_gene_dict.get)] == 1:
        if sec_essentiality_dict[min(sec_essentiality_dict, key=sec_essentiality_dict.get)] == 1:
            secesscoregenes += sec_list_of_genes
        elif sec_essentiality_dict[max(sec_essentiality_dict, key=sec_essentiality_dict.get)] == 0:
            secnesscoregenes += sec_list_of_genes
        else:
            secsesscoregenes += sec_list_of_genes
    else:
        if not sum([sec_gene_dict[item] - sec_essentiality_dict[item] for item in sec_essentiality_dict.keys()]):
            secessaccessorygenes += sec_list_of_genes
        elif sec_essentiality_dict[max(sec_essentiality_dict, key=sec_essentiality_dict.get)] == 0:
            secnessaccessorygenes += sec_list_of_genes
        else:
            secsessaccessorygenes += sec_list_of_genes

    if ke_gene_dict[min(ke_gene_dict, key=ke_gene_dict.get)] == 1:
        if ke_essentiality_dict[min(ke_essentiality_dict, key=ke_essentiality_dict.get)] == 1:
            keesscoregenes += ke_list_of_genes
        elif ke_essentiality_dict[max(ke_essentiality_dict, key=ke_essentiality_dict.get)] == 0:
            kenesscoregenes += ke_list_of_genes
        else:
            kesesscoregenes += ke_list_of_genes
    else:
        if not sum([ke_gene_dict[item] - ke_essentiality_dict[item] for item in ke_essentiality_dict.keys()]):
            keessaccessorygenes += ke_list_of_genes
        elif ke_essentiality_dict[max(ke_essentiality_dict, key=ke_essentiality_dict.get)] == 0:
            kenessaccessorygenes += ke_list_of_genes
        else:
            kesessaccessorygenes += ke_list_of_genes

esscoregenes = list(set(esscoregenes))
esscoregenes.sort()
sesscoregenes = list(set(sesscoregenes))
sesscoregenes.sort()
nesscoregenes = list(set(nesscoregenes))
nesscoregenes.sort()
essaccessorygenes = list(set(essaccessorygenes))
essaccessorygenes.sort()
sessaccessorygenes = list(set(sessaccessorygenes))
sessaccessorygenes.sort()
nessaccessorygenes = list(set(nessaccessorygenes))
nessaccessorygenes.sort()

salesscoregenes = list(set(salesscoregenes))
salesscoregenes.sort()
salsesscoregenes = list(set(salsesscoregenes))
salsesscoregenes.sort()
salnesscoregenes = list(set(salnesscoregenes))
salnesscoregenes.sort()
salessaccessorygenes = list(set(salessaccessorygenes))
salessaccessorygenes.sort()
salsessaccessorygenes = list(set(salsessaccessorygenes))
salsessaccessorygenes.sort()
salnessaccessorygenes = list(set(salnessaccessorygenes))
salnessaccessorygenes.sort()

ecoesscoregenes = list(set(ecoesscoregenes))
ecoesscoregenes.sort()
ecosesscoregenes = list(set(ecosesscoregenes))
ecosesscoregenes.sort()
econesscoregenes = list(set(econesscoregenes))
econesscoregenes.sort()
ecoessaccessorygenes = list(set(ecoessaccessorygenes))
ecoessaccessorygenes.sort()
ecosessaccessorygenes = list(set(ecosessaccessorygenes))
ecosessaccessorygenes.sort()
econessaccessorygenes = list(set(econessaccessorygenes))
econessaccessorygenes.sort()

kleesscoregenes = list(set(kleesscoregenes))
kleesscoregenes.sort()
klesesscoregenes = list(set(klesesscoregenes))
klesesscoregenes.sort()
klenesscoregenes = list(set(klenesscoregenes))
klenesscoregenes.sort()
kleessaccessorygenes = list(set(kleessaccessorygenes))
kleessaccessorygenes.sort()
klesessaccessorygenes = list(set(klesessaccessorygenes))
klesessaccessorygenes.sort()
klenessaccessorygenes = list(set(klenessaccessorygenes))
klenessaccessorygenes.sort()

citesscoregenes = list(set(citesscoregenes))
citesscoregenes.sort()
citsesscoregenes = list(set(citsesscoregenes))
citsesscoregenes.sort()
citnesscoregenes = list(set(citnesscoregenes))
citnesscoregenes.sort()
citessaccessorygenes = list(set(citessaccessorygenes))
citessaccessorygenes.sort()
citsessaccessorygenes = list(set(citsessaccessorygenes))
citsessaccessorygenes.sort()
citnessaccessorygenes = list(set(citnessaccessorygenes))
citnessaccessorygenes.sort()

entesscoregenes = list(set(entesscoregenes))
entesscoregenes.sort()
entsesscoregenes = list(set(entsesscoregenes))
entsesscoregenes.sort()
entnesscoregenes = list(set(entnesscoregenes))
entnesscoregenes.sort()
entessaccessorygenes = list(set(entessaccessorygenes))
entessaccessorygenes.sort()
entsessaccessorygenes = list(set(entsessaccessorygenes))
entsessaccessorygenes.sort()
entnessaccessorygenes = list(set(entnessaccessorygenes))
entnessaccessorygenes.sort()

seesscoregenes = list(set(seesscoregenes))
seesscoregenes.sort()
sesesscoregenes = list(set(sesesscoregenes))
sesesscoregenes.sort()
senesscoregenes = list(set(senesscoregenes))
senesscoregenes.sort()
seessaccessorygenes = list(set(seessaccessorygenes))
seessaccessorygenes.sort()
sesessaccessorygenes = list(set(sesessaccessorygenes))
sesessaccessorygenes.sort()
senessaccessorygenes = list(set(senessaccessorygenes))
senessaccessorygenes.sort()

secesscoregenes = list(set(secesscoregenes))
secesscoregenes.sort()
secsesscoregenes = list(set(secsesscoregenes))
secsesscoregenes.sort()
secnesscoregenes = list(set(secnesscoregenes))
secnesscoregenes.sort()
secessaccessorygenes = list(set(secessaccessorygenes))
secessaccessorygenes.sort()
secsessaccessorygenes = list(set(secsessaccessorygenes))
secsessaccessorygenes.sort()
secnessaccessorygenes = list(set(secnessaccessorygenes))
secnessaccessorygenes.sort()

keesscoregenes = list(set(keesscoregenes))
keesscoregenes.sort()
kesesscoregenes = list(set(kesesscoregenes))
kesesscoregenes.sort()
kenesscoregenes = list(set(kenesscoregenes))
kenesscoregenes.sort()
keessaccessorygenes = list(set(keessaccessorygenes))
keessaccessorygenes.sort()
kesessaccessorygenes = list(set(kesessaccessorygenes))
kesessaccessorygenes.sort()
kenessaccessorygenes = list(set(kenessaccessorygenes))
kenessaccessorygenes.sort()

sequences = read_fasta_sequences(seqdb)

prev_name = ''
for item in esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'esscoregenesfile' in locals():
            esscoregenesfile.close()
        esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
    esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sesscoregenesfile' in locals():
            sesscoregenesfile.close()
        sesscoregenesfile = open(sesscoredir + '/' + name + '.fasta', 'w')
    sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'nesscoregenesfile' in locals():
            nesscoregenesfile.close()
        nesscoregenesfile = open(nesscoredir + '/' + name + '.fasta', 'w')
    nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'essaccessorygenesfile' in locals():
            essaccessorygenesfile.close()
        essaccessorygenesfile = open(essaccessorydir + '/' + name + '.fasta', 'w')
    essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sessaccessorygenesfile' in locals():
            sessaccessorygenesfile.close()
        sessaccessorygenesfile = open(sessaccessorydir + '/' + name + '.fasta', 'w')
    sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'nessaccessorygenesfile' in locals():
            nessaccessorygenesfile.close()
        nessaccessorygenesfile = open(nessaccessorydir + '/' + name + '.fasta', 'w')
    nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salesscoregenesfile' in locals():
            salesscoregenesfile.close()
        salesscoregenesfile = open(salesscoredir + '/' + name + '.fasta', 'w')
    salesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salsesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salsesscoregenesfile' in locals():
            salsesscoregenesfile.close()
        salsesscoregenesfile = open(salsesscoredir + '/' + name + '.fasta', 'w')
    salsesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salnesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salnesscoregenesfile' in locals():
            salnesscoregenesfile.close()
        salnesscoregenesfile = open(salnesscoredir + '/' + name + '.fasta', 'w')
    salnesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salessaccessorygenesfile' in locals():
            salessaccessorygenesfile.close()
        salessaccessorygenesfile = open(salessaccessorydir + '/' + name + '.fasta', 'w')
    salessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salsessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salsessaccessorygenesfile' in locals():
            salsessaccessorygenesfile.close()
        salsessaccessorygenesfile = open(salsessaccessorydir + '/' + name + '.fasta', 'w')
    salsessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in salnessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'salnessaccessorygenesfile' in locals():
            salnessaccessorygenesfile.close()
        salnessaccessorygenesfile = open(salnessaccessorydir + '/' + name + '.fasta', 'w')
    salnessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ecoesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ecoesscoregenesfile' in locals():
            ecoesscoregenesfile.close()
        ecoesscoregenesfile = open(ecoesscoredir + '/' + name + '.fasta', 'w')
    ecoesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ecosesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ecosesscoregenesfile' in locals():
            ecosesscoregenesfile.close()
        ecosesscoregenesfile = open(ecosesscoredir + '/' + name + '.fasta', 'w')
    ecosesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in econesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'econesscoregenesfile' in locals():
            econesscoregenesfile.close()
        econesscoregenesfile = open(econesscoredir + '/' + name + '.fasta', 'w')
    econesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ecoessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ecoessaccessorygenesfile' in locals():
            ecoessaccessorygenesfile.close()
        ecoessaccessorygenesfile = open(ecoessaccessorydir + '/' + name + '.fasta', 'w')
    ecoessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ecosessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ecosessaccessorygenesfile' in locals():
            ecosessaccessorygenesfile.close()
        ecosessaccessorygenesfile = open(ecosessaccessorydir + '/' + name + '.fasta', 'w')
    ecosessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in econessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'econessaccessorygenesfile' in locals():
            econessaccessorygenesfile.close()
        econessaccessorygenesfile = open(econessaccessorydir + '/' + name + '.fasta', 'w')
    econessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kleesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kleesscoregenesfile' in locals():
            kleesscoregenesfile.close()
        kleesscoregenesfile = open(kleesscoredir + '/' + name + '.fasta', 'w')
    kleesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in klesesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'klesesscoregenesfile' in locals():
            klesesscoregenesfile.close()
        klesesscoregenesfile = open(klesesscoredir + '/' + name + '.fasta', 'w')
    klesesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in klenesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'klenesscoregenesfile' in locals():
            klenesscoregenesfile.close()
        klenesscoregenesfile = open(klenesscoredir + '/' + name + '.fasta', 'w')
    klenesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kleessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kleessaccessorygenesfile' in locals():
            kleessaccessorygenesfile.close()
        kleessaccessorygenesfile = open(kleessaccessorydir + '/' + name + '.fasta', 'w')
    kleessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in klesessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'klesessaccessorygenesfile' in locals():
            klesessaccessorygenesfile.close()
        klesessaccessorygenesfile = open(klesessaccessorydir + '/' + name + '.fasta', 'w')
    klesessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in klenessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'klenessaccessorygenesfile' in locals():
            klenessaccessorygenesfile.close()
        klenessaccessorygenesfile = open(klenessaccessorydir + '/' + name + '.fasta', 'w')
    klenessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citesscoregenesfile' in locals():
            citesscoregenesfile.close()
        citesscoregenesfile = open(citesscoredir + '/' + name + '.fasta', 'w')
    citesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citsesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citsesscoregenesfile' in locals():
            citsesscoregenesfile.close()
        citsesscoregenesfile = open(citsesscoredir + '/' + name + '.fasta', 'w')
    citsesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citnesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citnesscoregenesfile' in locals():
            citnesscoregenesfile.close()
        citnesscoregenesfile = open(citnesscoredir + '/' + name + '.fasta', 'w')
    citnesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citessaccessorygenesfile' in locals():
            citessaccessorygenesfile.close()
        citessaccessorygenesfile = open(citessaccessorydir + '/' + name + '.fasta', 'w')
    citessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citsessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citsessaccessorygenesfile' in locals():
            citsessaccessorygenesfile.close()
        citsessaccessorygenesfile = open(citsessaccessorydir + '/' + name + '.fasta', 'w')
    citsessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in citnessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'citnessaccessorygenesfile' in locals():
            citnessaccessorygenesfile.close()
        citnessaccessorygenesfile = open(citnessaccessorydir + '/' + name + '.fasta', 'w')
    citnessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entesscoregenesfile' in locals():
            entesscoregenesfile.close()
        entesscoregenesfile = open(entesscoredir + '/' + name + '.fasta', 'w')
    entesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entsesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entsesscoregenesfile' in locals():
            entsesscoregenesfile.close()
        entsesscoregenesfile = open(entsesscoredir + '/' + name + '.fasta', 'w')
    entsesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entnesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entnesscoregenesfile' in locals():
            entnesscoregenesfile.close()
        entnesscoregenesfile = open(entnesscoredir + '/' + name + '.fasta', 'w')
    entnesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entessaccessorygenesfile' in locals():
            entessaccessorygenesfile.close()
        entessaccessorygenesfile = open(entessaccessorydir + '/' + name + '.fasta', 'w')
    entessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entsessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entsessaccessorygenesfile' in locals():
            entsessaccessorygenesfile.close()
        entsessaccessorygenesfile = open(entsessaccessorydir + '/' + name + '.fasta', 'w')
    entsessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in entnessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'entnessaccessorygenesfile' in locals():
            entnessaccessorygenesfile.close()
        entnessaccessorygenesfile = open(entnessaccessorydir + '/' + name + '.fasta', 'w')
    entnessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in seesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'seesscoregenesfile' in locals():
            seesscoregenesfile.close()
        seesscoregenesfile = open(seesscoredir + '/' + name + '.fasta', 'w')
    seesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sesesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sesesscoregenesfile' in locals():
            sesesscoregenesfile.close()
        sesesscoregenesfile = open(sesesscoredir + '/' + name + '.fasta', 'w')
    sesesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in senesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'senesscoregenesfile' in locals():
            senesscoregenesfile.close()
        senesscoregenesfile = open(senesscoredir + '/' + name + '.fasta', 'w')
    senesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in seessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'seessaccessorygenesfile' in locals():
            seessaccessorygenesfile.close()
        seessaccessorygenesfile = open(seessaccessorydir + '/' + name + '.fasta', 'w')
    seessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sesessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sesessaccessorygenesfile' in locals():
            sesessaccessorygenesfile.close()
        sesessaccessorygenesfile = open(sesessaccessorydir + '/' + name + '.fasta', 'w')
    sesessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in senessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'senessaccessorygenesfile' in locals():
            senessaccessorygenesfile.close()
        senessaccessorygenesfile = open(senessaccessorydir + '/' + name + '.fasta', 'w')
    senessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secesscoregenesfile' in locals():
            secesscoregenesfile.close()
        secesscoregenesfile = open(secesscoredir + '/' + name + '.fasta', 'w')
    secesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secsesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secsesscoregenesfile' in locals():
            secsesscoregenesfile.close()
        secsesscoregenesfile = open(secsesscoredir + '/' + name + '.fasta', 'w')
    secsesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secnesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secnesscoregenesfile' in locals():
            secnesscoregenesfile.close()
        secnesscoregenesfile = open(secnesscoredir + '/' + name + '.fasta', 'w')
    secnesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secessaccessorygenesfile' in locals():
            secessaccessorygenesfile.close()
        secessaccessorygenesfile = open(secessaccessorydir + '/' + name + '.fasta', 'w')
    secessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secsessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secsessaccessorygenesfile' in locals():
            secsessaccessorygenesfile.close()
        secsessaccessorygenesfile = open(secsessaccessorydir + '/' + name + '.fasta', 'w')
    secsessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in secnessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'secnessaccessorygenesfile' in locals():
            secnessaccessorygenesfile.close()
        secnessaccessorygenesfile = open(secnessaccessorydir + '/' + name + '.fasta', 'w')
    secnessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in keesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'keesscoregenesfile' in locals():
            keesscoregenesfile.close()
        keesscoregenesfile = open(keesscoredir + '/' + name + '.fasta', 'w')
    keesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kesesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kesesscoregenesfile' in locals():
            kesesscoregenesfile.close()
        kesesscoregenesfile = open(kesesscoredir + '/' + name + '.fasta', 'w')
    kesesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kenesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kenesscoregenesfile' in locals():
            kenesscoregenesfile.close()
        kenesscoregenesfile = open(kenesscoredir + '/' + name + '.fasta', 'w')
    kenesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in keessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'keessaccessorygenesfile' in locals():
            keessaccessorygenesfile.close()
        keessaccessorygenesfile = open(keessaccessorydir + '/' + name + '.fasta', 'w')
    keessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kesessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kesessaccessorygenesfile' in locals():
            kesessaccessorygenesfile.close()
        kesessaccessorygenesfile = open(kesessaccessorydir + '/' + name + '.fasta', 'w')
    kesessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in kenessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'kenessaccessorygenesfile' in locals():
            kenessaccessorygenesfile.close()
        kenessaccessorygenesfile = open(kenessaccessorydir + '/' + name + '.fasta', 'w')
    kenessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name