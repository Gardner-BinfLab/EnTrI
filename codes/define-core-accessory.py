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

ty2esscoredir = outdir + '/salmonellatyphi-core-essential-genomes'
ty2sesscoredir = outdir + '/salmonellatyphi-core-sometimes-essential-genomes'
ty2nesscoredir = outdir + '/salmonellatyphi-core-never-essential-genomes'
ty2essaccessorydir = outdir + '/salmonellatyphi-accessory-essential-genomes'
ty2sessaccessorydir = outdir + '/salmonellatyphi-accessory-sometimes-essential-genomes'
ty2nessaccessorydir = outdir + '/salmonellatyphi-accessory-never-essential-genomes'

p12esscoredir = outdir + '/salmonellap125109-core-essential-genomes'
p12sesscoredir = outdir + '/salmonellap125109-core-sometimes-essential-genomes'
p12nesscoredir = outdir + '/salmonellap125109-core-never-essential-genomes'
p12essaccessorydir = outdir + '/salmonellap125109-accessory-essential-genomes'
p12sessaccessorydir = outdir + '/salmonellap125109-accessory-sometimes-essential-genomes'
p12nessaccessorydir = outdir + '/salmonellap125109-accessory-never-essential-genomes'

sl1esscoredir = outdir + '/salmonellasl1344-core-essential-genomes'
sl1sesscoredir = outdir + '/salmonellasl1344-core-sometimes-essential-genomes'
sl1nesscoredir = outdir + '/salmonellasl1344-core-never-essential-genomes'
sl1essaccessorydir = outdir + '/salmonellasl1344-accessory-essential-genomes'
sl1sessaccessorydir = outdir + '/salmonellasl1344-accessory-sometimes-essential-genomes'
sl1nessaccessorydir = outdir + '/salmonellasl1344-accessory-never-essential-genomes'

a13esscoredir = outdir + '/salmonellaa130-core-essential-genomes'
a13sesscoredir = outdir + '/salmonellaa130-core-sometimes-essential-genomes'
a13nesscoredir = outdir + '/salmonellaa130-core-never-essential-genomes'
a13essaccessorydir = outdir + '/salmonellaa130-accessory-essential-genomes'
a13sessaccessorydir = outdir + '/salmonellaa130-accessory-sometimes-essential-genomes'
a13nessaccessorydir = outdir + '/salmonellaa130-accessory-never-essential-genomes'

d23esscoredir = outdir + '/salmonellad23580-core-essential-genomes'
d23sesscoredir = outdir + '/salmonellad23580-core-sometimes-essential-genomes'
d23nesscoredir = outdir + '/salmonellad23580-core-never-essential-genomes'
d23essaccessorydir = outdir + '/salmonellad23580-accessory-essential-genomes'
d23sessaccessorydir = outdir + '/salmonellad23580-accessory-sometimes-essential-genomes'
d23nessaccessorydir = outdir + '/salmonellad23580-accessory-never-essential-genomes'

st1esscoredir = outdir + '/ecolist131-core-essential-genomes'
st1sesscoredir = outdir + '/ecolist131-core-sometimes-essential-genomes'
st1nesscoredir = outdir + '/ecolist131-core-never-essential-genomes'
st1essaccessorydir = outdir + '/ecolist131-accessory-essential-genomes'
st1sessaccessorydir = outdir + '/ecolist131-accessory-sometimes-essential-genomes'
st1nessaccessorydir = outdir + '/ecolist131-accessory-never-essential-genomes'

cs1esscoredir = outdir + '/ecolics17-core-essential-genomes'
cs1sesscoredir = outdir + '/ecolics17-core-sometimes-essential-genomes'
cs1nesscoredir = outdir + '/scolics17-core-never-essential-genomes'
cs1essaccessorydir = outdir + '/ecolics17-accessory-essential-genomes'
cs1sessaccessorydir = outdir + '/ecolics17-accessory-sometimes-essential-genomes'
cs1nessaccessorydir = outdir + '/ecolics17-accessory-never-essential-genomes'

h10esscoredir = outdir + '/ecolih10407-core-essential-genomes'
h10sesscoredir = outdir + '/ecolih10407-core-sometimes-essential-genomes'
h10nesscoredir = outdir + '/ecolih10407-core-never-essential-genomes'
h10essaccessorydir = outdir + '/ecolih10407-accessory-essential-genomes'
h10sessaccessorydir = outdir + '/ecolih10407-accessory-sometimes-essential-genomes'
h10nessaccessorydir = outdir + '/ecolih10407-accessory-never-essential-genomes'

rh2esscoredir = outdir + '/klebsiellarh201207-core-essential-genomes'
rh2sesscoredir = outdir + '/klebsiellarh201207-core-sometimes-essential-genomes'
rh2nesscoredir = outdir + '/klebsiellarh201207-core-never-essential-genomes'
rh2essaccessorydir = outdir + '/klebsiellarh201207-accessory-essential-genomes'
rh2sessaccessorydir = outdir + '/klebsiellarh201207-accessory-sometimes-essential-genomes'
rh2nessaccessorydir = outdir + '/klebsiellarh201207-accessory-never-essential-genomes'

eclesscoredir = outdir + '/klebsiellaecl8-core-essential-genomes'
eclsesscoredir = outdir + '/klebsiellaecl8-core-sometimes-essential-genomes'
eclnesscoredir = outdir + '/klebsiellaecl8-core-never-essential-genomes'
eclessaccessorydir = outdir + '/klebsiellaecl8-accessory-essential-genomes'
eclsessaccessorydir = outdir + '/klebsiellaecl8-accessory-sometimes-essential-genomes'
eclnessaccessorydir = outdir + '/klebsiellaecl8-accessory-never-essential-genomes'

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

makedir(ty2esscoredir)
makedir(ty2sesscoredir)
makedir(ty2nesscoredir)
makedir(ty2essaccessorydir)
makedir(ty2sessaccessorydir)
makedir(ty2nessaccessorydir)

makedir(p12esscoredir)
makedir(p12sesscoredir)
makedir(p12nesscoredir)
makedir(p12essaccessorydir)
makedir(p12sessaccessorydir)
makedir(p12nessaccessorydir)

makedir(sl1esscoredir)
makedir(sl1sesscoredir)
makedir(sl1nesscoredir)
makedir(sl1essaccessorydir)
makedir(sl1sessaccessorydir)
makedir(sl1nessaccessorydir)

makedir(a13esscoredir)
makedir(a13sesscoredir)
makedir(a13nesscoredir)
makedir(a13essaccessorydir)
makedir(a13sessaccessorydir)
makedir(a13nessaccessorydir)

makedir(d23esscoredir)
makedir(d23sesscoredir)
makedir(d23nesscoredir)
makedir(d23essaccessorydir)
makedir(d23sessaccessorydir)
makedir(d23nessaccessorydir)

makedir(st1esscoredir)
makedir(st1sesscoredir)
makedir(st1nesscoredir)
makedir(st1essaccessorydir)
makedir(st1sessaccessorydir)
makedir(st1nessaccessorydir)

makedir(cs1esscoredir)
makedir(cs1sesscoredir)
makedir(cs1nesscoredir)
makedir(cs1essaccessorydir)
makedir(cs1sessaccessorydir)
makedir(cs1nessaccessorydir)

makedir(h10esscoredir)
makedir(h10sesscoredir)
makedir(h10nesscoredir)
makedir(h10essaccessorydir)
makedir(h10sessaccessorydir)
makedir(h10nessaccessorydir)

makedir(rh2esscoredir)
makedir(rh2sesscoredir)
makedir(rh2nesscoredir)
makedir(rh2essaccessorydir)
makedir(rh2sessaccessorydir)
makedir(rh2nessaccessorydir)

makedir(eclesscoredir)
makedir(eclsesscoredir)
makedir(eclnesscoredir)
makedir(eclessaccessorydir)
makedir(eclsessaccessorydir)
makedir(eclnessaccessorydir)

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

ty2esscoregenes = []
ty2sesscoregenes = []
ty2nesscoregenes = []
ty2essaccessorygenes = []
ty2sessaccessorygenes = []
ty2nessaccessorygenes = []

p12esscoregenes = []
p12sesscoregenes = []
p12nesscoregenes = []
p12essaccessorygenes = []
p12sessaccessorygenes = []
p12nessaccessorygenes = []

sl1esscoregenes = []
sl1sesscoregenes = []
sl1nesscoregenes = []
sl1essaccessorygenes = []
sl1sessaccessorygenes = []
sl1nessaccessorygenes = []

a13esscoregenes = []
a13sesscoregenes = []
a13nesscoregenes = []
a13essaccessorygenes = []
a13sessaccessorygenes = []
a13nessaccessorygenes = []

d23esscoregenes = []
d23sesscoregenes = []
d23nesscoregenes = []
d23essaccessorygenes = []
d23sessaccessorygenes = []
d23nessaccessorygenes = []

st1esscoregenes = []
st1sesscoregenes = []
st1nesscoregenes = []
st1essaccessorygenes = []
st1sessaccessorygenes = []
st1nessaccessorygenes = []

cs1esscoregenes = []
cs1sesscoregenes = []
cs1nesscoregenes = []
cs1essaccessorygenes = []
cs1sessaccessorygenes = []
cs1nessaccessorygenes = []

h10esscoregenes = []
h10sesscoregenes = []
h10nesscoregenes = []
h10essaccessorygenes = []
h10sessaccessorygenes = []
h10nessaccessorygenes = []

rh2esscoregenes = []
rh2sesscoregenes = []
rh2nesscoregenes = []
rh2essaccessorygenes = []
rh2sessaccessorygenes = []
rh2nessaccessorygenes = []

eclesscoregenes = []
eclsesscoregenes = []
eclnesscoregenes = []
eclessaccessorygenes = []
eclsessaccessorygenes = []
eclnessaccessorygenes = []

species_names = ["BN373", "CS17", "ENC", "ERS227112", "ETEC", "NCTC13441", "ROD", "SEN", "SL1344", "STM", "STMMW", "t"]
num_species = len(species_names)
gene_dict = {species_names[i]: 0 for i in range(num_species)}
essentiality_dict = {species_names[i]: 0 for i in range(num_species)}

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

ty2_species_names = ["t"]
ty2_num_species = len(ty2_species_names)
ty2_gene_dict = {ty2_species_names[i]: 0 for i in range(ty2_num_species)}
ty2_essentiality_dict = {ty2_species_names[i]: 0 for i in range(ty2_num_species)}

p12_species_names = ["SEN"]
p12_num_species = len(p12_species_names)
p12_gene_dict = {p12_species_names[i]: 0 for i in range(p12_num_species)}
p12_essentiality_dict = {p12_species_names[i]: 0 for i in range(p12_num_species)}

sl1_species_names = ["SL1344"]
sl1_num_species = len(sl1_species_names)
sl1_gene_dict = {sl1_species_names[i]: 0 for i in range(sl1_num_species)}
sl1_essentiality_dict = {sl1_species_names[i]: 0 for i in range(sl1_num_species)}

a13_species_names = ["STM"]
a13_num_species = len(a13_species_names)
a13_gene_dict = {a13_species_names[i]: 0 for i in range(a13_num_species)}
a13_essentiality_dict = {a13_species_names[i]: 0 for i in range(a13_num_species)}

d23_species_names = ["STMMW"]
d23_num_species = len(d23_species_names)
d23_gene_dict = {d23_species_names[i]: 0 for i in range(d23_num_species)}
d23_essentiality_dict = {d23_species_names[i]: 0 for i in range(d23_num_species)}

st1_species_names = ["NCTC13441"]
st1_num_species = len(st1_species_names)
st1_gene_dict = {st1_species_names[i]: 0 for i in range(st1_num_species)}
st1_essentiality_dict = {st1_species_names[i]: 0 for i in range(st1_num_species)}

cs1_species_names = ["CS17"]
cs1_num_species = len(cs1_species_names)
cs1_gene_dict = {cs1_species_names[i]: 0 for i in range(cs1_num_species)}
cs1_essentiality_dict = {cs1_species_names[i]: 0 for i in range(cs1_num_species)}

h10_species_names = ["ETEC"]
h10_num_species = len(h10_species_names)
h10_gene_dict = {h10_species_names[i]: 0 for i in range(h10_num_species)}
h10_essentiality_dict = {h10_species_names[i]: 0 for i in range(h10_num_species)}

rh2_species_names = ["ERS227112"]
rh2_num_species = len(rh2_species_names)
rh2_gene_dict = {rh2_species_names[i]: 0 for i in range(rh2_num_species)}
rh2_essentiality_dict = {rh2_species_names[i]: 0 for i in range(rh2_num_species)}

ecl_species_names = ["BN373"]
ecl_num_species = len(ecl_species_names)
ecl_gene_dict = {ecl_species_names[i]: 0 for i in range(ecl_num_species)}
ecl_essentiality_dict = {ecl_species_names[i]: 0 for i in range(ecl_num_species)}

list_of_files = listdir(clusters)
for filename in list_of_files:
    for key in gene_dict.keys():
        gene_dict[key] = 0
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

    for key in ty2_gene_dict.keys():
        ty2_gene_dict[key] = 0
        ty2_essentiality_dict[key] = 0

    for key in p12_gene_dict.keys():
        p12_gene_dict[key] = 0
        p12_essentiality_dict[key] = 0

    for key in sl1_gene_dict.keys():
        sl1_gene_dict[key] = 0
        sl1_essentiality_dict[key] = 0

    for key in a13_gene_dict.keys():
        a13_gene_dict[key] = 0
        a13_essentiality_dict[key] = 0

    for key in d23_gene_dict.keys():
        d23_gene_dict[key] = 0
        d23_essentiality_dict[key] = 0

    for key in st1_gene_dict.keys():
        st1_gene_dict[key] = 0
        st1_essentiality_dict[key] = 0

    for key in cs1_gene_dict.keys():
        cs1_gene_dict[key] = 0
        cs1_essentiality_dict[key] = 0

    for key in h10_gene_dict.keys():
        h10_gene_dict[key] = 0
        h10_essentiality_dict[key] = 0

    for key in rh2_gene_dict.keys():
        rh2_gene_dict[key] = 0
        rh2_essentiality_dict[key] = 0

    for key in ecl_gene_dict.keys():
        ecl_gene_dict[key] = 0
        ecl_essentiality_dict[key] = 0

    list_of_genes = []
    sal_list_of_genes = []
    eco_list_of_genes = []
    kle_list_of_genes = []
    cit_list_of_genes = []
    ent_list_of_genes = []
    se_list_of_genes = []
    sec_list_of_genes = []
    ke_list_of_genes = []
    ty2_list_of_genes = []
    p12_list_of_genes = []
    sl1_list_of_genes = []
    a13_list_of_genes = []
    d23_list_of_genes = []
    st1_list_of_genes = []
    cs1_list_of_genes = []
    h10_list_of_genes = []
    rh2_list_of_genes = []
    ecl_list_of_genes = []

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
                    if float(cells[4]) < 0.2:
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

                if name in ty2_gene_dict.keys():
                    ty2_list_of_genes.append(cells[1])
                    ty2_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        ty2_essentiality_dict[name] = 1

                if name in p12_gene_dict.keys():
                    p12_list_of_genes.append(cells[1])
                    p12_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        p12_essentiality_dict[name] = 1

                if name in sl1_gene_dict.keys():
                    sl1_list_of_genes.append(cells[1])
                    sl1_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        sl1_essentiality_dict[name] = 1

                if name in a13_gene_dict.keys():
                    a13_list_of_genes.append(cells[1])
                    a13_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        a13_essentiality_dict[name] = 1

                if name in d23_gene_dict.keys():
                    d23_list_of_genes.append(cells[1])
                    d23_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        d23_essentiality_dict[name] = 1

                if name in st1_gene_dict.keys():
                    st1_list_of_genes.append(cells[1])
                    st1_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        st1_essentiality_dict[name] = 1

                if name in cs1_gene_dict.keys():
                    cs1_list_of_genes.append(cells[1])
                    cs1_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        cs1_essentiality_dict[name] = 1

                if name in h10_gene_dict.keys():
                    h10_list_of_genes.append(cells[1])
                    h10_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        h10_essentiality_dict[name] = 1

                if name in rh2_gene_dict.keys():
                    rh2_list_of_genes.append(cells[1])
                    rh2_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        rh2_essentiality_dict[name] = 1

                if name in ecl_gene_dict.keys():
                    ecl_list_of_genes.append(cells[1])
                    ecl_gene_dict[name] = 1
                    if float(cells[4]) < 0.2:
                        ecl_essentiality_dict[name] = 1

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

    if ty2_gene_dict[min(ty2_gene_dict, key=ty2_gene_dict.get)] == 1:
        if ty2_essentiality_dict[min(ty2_essentiality_dict, key=ty2_essentiality_dict.get)] == 1:
            ty2esscoregenes += ty2_list_of_genes
        elif ty2_essentiality_dict[max(ty2_essentiality_dict, key=ty2_essentiality_dict.get)] == 0:
            ty2nesscoregenes += ty2_list_of_genes
        else:
            ty2sesscoregenes += ty2_list_of_genes
    else:
        if not sum([ty2_gene_dict[item] - ty2_essentiality_dict[item] for item in ty2_essentiality_dict.keys()]):
            ty2essaccessorygenes += ty2_list_of_genes
        elif ty2_essentiality_dict[max(ty2_essentiality_dict, key=ty2_essentiality_dict.get)] == 0:
            ty2nessaccessorygenes += ty2_list_of_genes
        else:
            ty2sessaccessorygenes += ty2_list_of_genes

    if p12_gene_dict[min(p12_gene_dict, key=p12_gene_dict.get)] == 1:
        if p12_essentiality_dict[min(p12_essentiality_dict, key=p12_essentiality_dict.get)] == 1:
            p12esscoregenes += p12_list_of_genes
        elif p12_essentiality_dict[max(p12_essentiality_dict, key=p12_essentiality_dict.get)] == 0:
            p12nesscoregenes += p12_list_of_genes
        else:
            p12sesscoregenes += p12_list_of_genes
    else:
        if not sum([p12_gene_dict[item] - p12_essentiality_dict[item] for item in p12_essentiality_dict.keys()]):
            p12essaccessorygenes += p12_list_of_genes
        elif p12_essentiality_dict[max(p12_essentiality_dict, key=p12_essentiality_dict.get)] == 0:
            p12nessaccessorygenes += p12_list_of_genes
        else:
            p12sessaccessorygenes += p12_list_of_genes

    if sl1_gene_dict[min(sl1_gene_dict, key=sl1_gene_dict.get)] == 1:
        if sl1_essentiality_dict[min(sl1_essentiality_dict, key=sl1_essentiality_dict.get)] == 1:
            sl1esscoregenes += sl1_list_of_genes
        elif sl1_essentiality_dict[max(sl1_essentiality_dict, key=sl1_essentiality_dict.get)] == 0:
            sl1nesscoregenes += sl1_list_of_genes
        else:
            sl1sesscoregenes += sl1_list_of_genes
    else:
        if not sum([sl1_gene_dict[item] - sl1_essentiality_dict[item] for item in sl1_essentiality_dict.keys()]):
            sl1essaccessorygenes += sl1_list_of_genes
        elif sl1_essentiality_dict[max(sl1_essentiality_dict, key=sl1_essentiality_dict.get)] == 0:
            sl1nessaccessorygenes += sl1_list_of_genes
        else:
            sl1sessaccessorygenes += sl1_list_of_genes

    if a13_gene_dict[min(a13_gene_dict, key=a13_gene_dict.get)] == 1:
        if a13_essentiality_dict[min(a13_essentiality_dict, key=a13_essentiality_dict.get)] == 1:
            a13esscoregenes += a13_list_of_genes
        elif a13_essentiality_dict[max(a13_essentiality_dict, key=a13_essentiality_dict.get)] == 0:
            a13nesscoregenes += a13_list_of_genes
        else:
            a13sesscoregenes += a13_list_of_genes
    else:
        if not sum([a13_gene_dict[item] - a13_essentiality_dict[item] for item in a13_essentiality_dict.keys()]):
            a13essaccessorygenes += a13_list_of_genes
        elif a13_essentiality_dict[max(a13_essentiality_dict, key=a13_essentiality_dict.get)] == 0:
            a13nessaccessorygenes += a13_list_of_genes
        else:
            a13sessaccessorygenes += a13_list_of_genes

    if d23_gene_dict[min(d23_gene_dict, key=d23_gene_dict.get)] == 1:
        if d23_essentiality_dict[min(d23_essentiality_dict, key=d23_essentiality_dict.get)] == 1:
            d23esscoregenes += d23_list_of_genes
        elif d23_essentiality_dict[max(d23_essentiality_dict, key=d23_essentiality_dict.get)] == 0:
            d23nesscoregenes += d23_list_of_genes
        else:
            d23sesscoregenes += d23_list_of_genes
    else:
        if not sum([d23_gene_dict[item] - d23_essentiality_dict[item] for item in d23_essentiality_dict.keys()]):
            d23essaccessorygenes += d23_list_of_genes
        elif d23_essentiality_dict[max(d23_essentiality_dict, key=d23_essentiality_dict.get)] == 0:
            d23nessaccessorygenes += d23_list_of_genes
        else:
            d23sessaccessorygenes += d23_list_of_genes

    if st1_gene_dict[min(st1_gene_dict, key=st1_gene_dict.get)] == 1:
        if st1_essentiality_dict[min(st1_essentiality_dict, key=st1_essentiality_dict.get)] == 1:
            st1esscoregenes += st1_list_of_genes
        elif st1_essentiality_dict[max(st1_essentiality_dict, key=st1_essentiality_dict.get)] == 0:
            st1nesscoregenes += st1_list_of_genes
        else:
            st1sesscoregenes += st1_list_of_genes
    else:
        if not sum([st1_gene_dict[item] - st1_essentiality_dict[item] for item in st1_essentiality_dict.keys()]):
            st1essaccessorygenes += st1_list_of_genes
        elif st1_essentiality_dict[max(st1_essentiality_dict, key=st1_essentiality_dict.get)] == 0:
            st1nessaccessorygenes += st1_list_of_genes
        else:
            st1sessaccessorygenes += st1_list_of_genes

    if cs1_gene_dict[min(cs1_gene_dict, key=cs1_gene_dict.get)] == 1:
        if cs1_essentiality_dict[min(cs1_essentiality_dict, key=cs1_essentiality_dict.get)] == 1:
            cs1esscoregenes += cs1_list_of_genes
        elif cs1_essentiality_dict[max(cs1_essentiality_dict, key=cs1_essentiality_dict.get)] == 0:
            cs1nesscoregenes += cs1_list_of_genes
        else:
            cs1sesscoregenes += cs1_list_of_genes
    else:
        if not sum([cs1_gene_dict[item] - cs1_essentiality_dict[item] for item in cs1_essentiality_dict.keys()]):
            cs1essaccessorygenes += cs1_list_of_genes
        elif cs1_essentiality_dict[max(cs1_essentiality_dict, key=cs1_essentiality_dict.get)] == 0:
            cs1nessaccessorygenes += cs1_list_of_genes
        else:
            cs1sessaccessorygenes += cs1_list_of_genes

    if h10_gene_dict[min(h10_gene_dict, key=h10_gene_dict.get)] == 1:
        if h10_essentiality_dict[min(h10_essentiality_dict, key=h10_essentiality_dict.get)] == 1:
            h10esscoregenes += h10_list_of_genes
        elif h10_essentiality_dict[max(h10_essentiality_dict, key=h10_essentiality_dict.get)] == 0:
            h10nesscoregenes += h10_list_of_genes
        else:
            h10sesscoregenes += h10_list_of_genes
    else:
        if not sum([h10_gene_dict[item] - h10_essentiality_dict[item] for item in h10_essentiality_dict.keys()]):
            h10essaccessorygenes += h10_list_of_genes
        elif h10_essentiality_dict[max(h10_essentiality_dict, key=h10_essentiality_dict.get)] == 0:
            h10nessaccessorygenes += h10_list_of_genes
        else:
            h10sessaccessorygenes += h10_list_of_genes

    if rh2_gene_dict[min(rh2_gene_dict, key=rh2_gene_dict.get)] == 1:
        if rh2_essentiality_dict[min(rh2_essentiality_dict, key=rh2_essentiality_dict.get)] == 1:
            rh2esscoregenes += rh2_list_of_genes
        elif rh2_essentiality_dict[max(rh2_essentiality_dict, key=rh2_essentiality_dict.get)] == 0:
            rh2nesscoregenes += rh2_list_of_genes
        else:
            rh2sesscoregenes += rh2_list_of_genes
    else:
        if not sum([rh2_gene_dict[item] - rh2_essentiality_dict[item] for item in rh2_essentiality_dict.keys()]):
            rh2essaccessorygenes += rh2_list_of_genes
        elif rh2_essentiality_dict[max(rh2_essentiality_dict, key=rh2_essentiality_dict.get)] == 0:
            rh2nessaccessorygenes += rh2_list_of_genes
        else:
            rh2sessaccessorygenes += rh2_list_of_genes

    if ecl_gene_dict[min(ecl_gene_dict, key=ecl_gene_dict.get)] == 1:
        if ecl_essentiality_dict[min(ecl_essentiality_dict, key=ecl_essentiality_dict.get)] == 1:
            eclesscoregenes += ecl_list_of_genes
        elif ecl_essentiality_dict[max(ecl_essentiality_dict, key=ecl_essentiality_dict.get)] == 0:
            eclnesscoregenes += ecl_list_of_genes
        else:
            eclsesscoregenes += ecl_list_of_genes
    else:
        if not sum([ecl_gene_dict[item] - ecl_essentiality_dict[item] for item in ecl_essentiality_dict.keys()]):
            eclessaccessorygenes += ecl_list_of_genes
        elif ecl_essentiality_dict[max(ecl_essentiality_dict, key=ecl_essentiality_dict.get)] == 0:
            eclnessaccessorygenes += ecl_list_of_genes
        else:
            eclsessaccessorygenes += ecl_list_of_genes

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

ty2esscoregenes = list(set(ty2esscoregenes))
ty2esscoregenes.sort()
ty2sesscoregenes = list(set(ty2sesscoregenes))
ty2sesscoregenes.sort()
ty2nesscoregenes = list(set(ty2nesscoregenes))
ty2nesscoregenes.sort()
ty2essaccessorygenes = list(set(ty2essaccessorygenes))
ty2essaccessorygenes.sort()
ty2sessaccessorygenes = list(set(ty2sessaccessorygenes))
ty2sessaccessorygenes.sort()
ty2nessaccessorygenes = list(set(ty2nessaccessorygenes))
ty2nessaccessorygenes.sort()

p12esscoregenes = list(set(p12esscoregenes))
p12esscoregenes.sort()
p12sesscoregenes = list(set(p12sesscoregenes))
p12sesscoregenes.sort()
p12nesscoregenes = list(set(p12nesscoregenes))
p12nesscoregenes.sort()
p12essaccessorygenes = list(set(p12essaccessorygenes))
p12essaccessorygenes.sort()
p12sessaccessorygenes = list(set(p12sessaccessorygenes))
p12sessaccessorygenes.sort()
p12nessaccessorygenes = list(set(p12nessaccessorygenes))
p12nessaccessorygenes.sort()

sl1esscoregenes = list(set(sl1esscoregenes))
sl1esscoregenes.sort()
sl1sesscoregenes = list(set(sl1sesscoregenes))
sl1sesscoregenes.sort()
sl1nesscoregenes = list(set(sl1nesscoregenes))
sl1nesscoregenes.sort()
sl1essaccessorygenes = list(set(sl1essaccessorygenes))
sl1essaccessorygenes.sort()
sl1sessaccessorygenes = list(set(sl1sessaccessorygenes))
sl1sessaccessorygenes.sort()
sl1nessaccessorygenes = list(set(sl1nessaccessorygenes))
sl1nessaccessorygenes.sort()

a13esscoregenes = list(set(a13esscoregenes))
a13esscoregenes.sort()
a13sesscoregenes = list(set(a13sesscoregenes))
a13sesscoregenes.sort()
a13nesscoregenes = list(set(a13nesscoregenes))
a13nesscoregenes.sort()
a13essaccessorygenes = list(set(a13essaccessorygenes))
a13essaccessorygenes.sort()
a13sessaccessorygenes = list(set(a13sessaccessorygenes))
a13sessaccessorygenes.sort()
a13nessaccessorygenes = list(set(a13nessaccessorygenes))
a13nessaccessorygenes.sort()

d23esscoregenes = list(set(d23esscoregenes))
d23esscoregenes.sort()
d23sesscoregenes = list(set(d23sesscoregenes))
d23sesscoregenes.sort()
d23nesscoregenes = list(set(d23nesscoregenes))
d23nesscoregenes.sort()
d23essaccessorygenes = list(set(d23essaccessorygenes))
d23essaccessorygenes.sort()
d23sessaccessorygenes = list(set(d23sessaccessorygenes))
d23sessaccessorygenes.sort()
d23nessaccessorygenes = list(set(d23nessaccessorygenes))
d23nessaccessorygenes.sort()

st1esscoregenes = list(set(st1esscoregenes))
st1esscoregenes.sort()
st1sesscoregenes = list(set(st1sesscoregenes))
st1sesscoregenes.sort()
st1nesscoregenes = list(set(st1nesscoregenes))
st1nesscoregenes.sort()
st1essaccessorygenes = list(set(st1essaccessorygenes))
st1essaccessorygenes.sort()
st1sessaccessorygenes = list(set(st1sessaccessorygenes))
st1sessaccessorygenes.sort()
st1nessaccessorygenes = list(set(st1nessaccessorygenes))
st1nessaccessorygenes.sort()

cs1esscoregenes = list(set(cs1esscoregenes))
cs1esscoregenes.sort()
cs1sesscoregenes = list(set(cs1sesscoregenes))
cs1sesscoregenes.sort()
cs1nesscoregenes = list(set(cs1nesscoregenes))
cs1nesscoregenes.sort()
cs1essaccessorygenes = list(set(cs1essaccessorygenes))
cs1essaccessorygenes.sort()
cs1sessaccessorygenes = list(set(cs1sessaccessorygenes))
cs1sessaccessorygenes.sort()
cs1nessaccessorygenes = list(set(cs1nessaccessorygenes))
cs1nessaccessorygenes.sort()

h10esscoregenes = list(set(h10esscoregenes))
h10esscoregenes.sort()
h10sesscoregenes = list(set(h10sesscoregenes))
h10sesscoregenes.sort()
h10nesscoregenes = list(set(h10nesscoregenes))
h10nesscoregenes.sort()
h10essaccessorygenes = list(set(h10essaccessorygenes))
h10essaccessorygenes.sort()
h10sessaccessorygenes = list(set(h10sessaccessorygenes))
h10sessaccessorygenes.sort()
h10nessaccessorygenes = list(set(h10nessaccessorygenes))
h10nessaccessorygenes.sort()

rh2esscoregenes = list(set(rh2esscoregenes))
rh2esscoregenes.sort()
rh2sesscoregenes = list(set(rh2sesscoregenes))
rh2sesscoregenes.sort()
rh2nesscoregenes = list(set(rh2nesscoregenes))
rh2nesscoregenes.sort()
rh2essaccessorygenes = list(set(rh2essaccessorygenes))
rh2essaccessorygenes.sort()
rh2sessaccessorygenes = list(set(rh2sessaccessorygenes))
rh2sessaccessorygenes.sort()
rh2nessaccessorygenes = list(set(rh2nessaccessorygenes))
rh2nessaccessorygenes.sort()

eclesscoregenes = list(set(eclesscoregenes))
eclesscoregenes.sort()
eclsesscoregenes = list(set(eclsesscoregenes))
eclsesscoregenes.sort()
eclnesscoregenes = list(set(eclnesscoregenes))
eclnesscoregenes.sort()
eclessaccessorygenes = list(set(eclessaccessorygenes))
eclessaccessorygenes.sort()
eclsessaccessorygenes = list(set(eclsessaccessorygenes))
eclsessaccessorygenes.sort()
eclnessaccessorygenes = list(set(eclnessaccessorygenes))
eclnessaccessorygenes.sort()

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

prev_name = ''
for item in ty2esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2esscoregenesfile' in locals():
            ty2esscoregenesfile.close()
        ty2esscoregenesfile = open(ty2esscoredir + '/' + name + '.fasta', 'w')
    ty2esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ty2sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2sesscoregenesfile' in locals():
            ty2sesscoregenesfile.close()
        ty2sesscoregenesfile = open(ty2sesscoredir + '/' + name + '.fasta', 'w')
    ty2sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ty2nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2nesscoregenesfile' in locals():
            ty2nesscoregenesfile.close()
        ty2nesscoregenesfile = open(ty2nesscoredir + '/' + name + '.fasta', 'w')
    ty2nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ty2essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2essaccessorygenesfile' in locals():
            ty2essaccessorygenesfile.close()
        ty2essaccessorygenesfile = open(ty2essaccessorydir + '/' + name + '.fasta', 'w')
    ty2essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ty2sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2sessaccessorygenesfile' in locals():
            ty2sessaccessorygenesfile.close()
        ty2sessaccessorygenesfile = open(ty2sessaccessorydir + '/' + name + '.fasta', 'w')
    ty2sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in ty2nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'ty2nessaccessorygenesfile' in locals():
            ty2nessaccessorygenesfile.close()
        ty2nessaccessorygenesfile = open(ty2nessaccessorydir + '/' + name + '.fasta', 'w')
    ty2nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

    prev_name = ''
for item in p12esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12esscoregenesfile' in locals():
            p12esscoregenesfile.close()
        p12esscoregenesfile = open(p12esscoredir + '/' + name + '.fasta', 'w')
    p12esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in p12sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12sesscoregenesfile' in locals():
            p12sesscoregenesfile.close()
        p12sesscoregenesfile = open(p12sesscoredir + '/' + name + '.fasta', 'w')
    p12sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in p12nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12nesscoregenesfile' in locals():
            p12nesscoregenesfile.close()
        p12nesscoregenesfile = open(p12nesscoredir + '/' + name + '.fasta', 'w')
    p12nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in p12essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12essaccessorygenesfile' in locals():
            p12essaccessorygenesfile.close()
        p12essaccessorygenesfile = open(p12essaccessorydir + '/' + name + '.fasta', 'w')
    p12essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in p12sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12sessaccessorygenesfile' in locals():
            p12sessaccessorygenesfile.close()
        p12sessaccessorygenesfile = open(p12sessaccessorydir + '/' + name + '.fasta', 'w')
    p12sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in p12nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'p12nessaccessorygenesfile' in locals():
            p12nessaccessorygenesfile.close()
        p12nessaccessorygenesfile = open(p12nessaccessorydir + '/' + name + '.fasta', 'w')
    p12nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1esscoregenesfile' in locals():
            sl1esscoregenesfile.close()
        sl1esscoregenesfile = open(sl1esscoredir + '/' + name + '.fasta', 'w')
    sl1esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1sesscoregenesfile' in locals():
            sl1sesscoregenesfile.close()
        sl1sesscoregenesfile = open(sl1sesscoredir + '/' + name + '.fasta', 'w')
    sl1sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1nesscoregenesfile' in locals():
            sl1nesscoregenesfile.close()
        sl1nesscoregenesfile = open(sl1nesscoredir + '/' + name + '.fasta', 'w')
    sl1nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1essaccessorygenesfile' in locals():
            sl1essaccessorygenesfile.close()
        sl1essaccessorygenesfile = open(sl1essaccessorydir + '/' + name + '.fasta', 'w')
    sl1essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1sessaccessorygenesfile' in locals():
            sl1sessaccessorygenesfile.close()
        sl1sessaccessorygenesfile = open(sl1sessaccessorydir + '/' + name + '.fasta', 'w')
    sl1sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in sl1nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'sl1nessaccessorygenesfile' in locals():
            sl1nessaccessorygenesfile.close()
        sl1nessaccessorygenesfile = open(sl1nessaccessorydir + '/' + name + '.fasta', 'w')
    sl1nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13esscoregenesfile' in locals():
            a13esscoregenesfile.close()
        a13esscoregenesfile = open(a13esscoredir + '/' + name + '.fasta', 'w')
    a13esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13sesscoregenesfile' in locals():
            a13sesscoregenesfile.close()
        a13sesscoregenesfile = open(a13sesscoredir + '/' + name + '.fasta', 'w')
    a13sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13nesscoregenesfile' in locals():
            a13nesscoregenesfile.close()
        a13nesscoregenesfile = open(a13nesscoredir + '/' + name + '.fasta', 'w')
    a13nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13essaccessorygenesfile' in locals():
            a13essaccessorygenesfile.close()
        a13essaccessorygenesfile = open(a13essaccessorydir + '/' + name + '.fasta', 'w')
    a13essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13sessaccessorygenesfile' in locals():
            a13sessaccessorygenesfile.close()
        a13sessaccessorygenesfile = open(a13sessaccessorydir + '/' + name + '.fasta', 'w')
    a13sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in a13nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'a13nessaccessorygenesfile' in locals():
            a13nessaccessorygenesfile.close()
        a13nessaccessorygenesfile = open(a13nessaccessorydir + '/' + name + '.fasta', 'w')
    a13nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23esscoregenesfile' in locals():
            d23esscoregenesfile.close()
        d23esscoregenesfile = open(d23esscoredir + '/' + name + '.fasta', 'w')
    d23esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23sesscoregenesfile' in locals():
            d23sesscoregenesfile.close()
        d23sesscoregenesfile = open(d23sesscoredir + '/' + name + '.fasta', 'w')
    d23sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23nesscoregenesfile' in locals():
            d23nesscoregenesfile.close()
        d23nesscoregenesfile = open(d23nesscoredir + '/' + name + '.fasta', 'w')
    d23nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23essaccessorygenesfile' in locals():
            d23essaccessorygenesfile.close()
        d23essaccessorygenesfile = open(d23essaccessorydir + '/' + name + '.fasta', 'w')
    d23essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23sessaccessorygenesfile' in locals():
            d23sessaccessorygenesfile.close()
        d23sessaccessorygenesfile = open(d23sessaccessorydir + '/' + name + '.fasta', 'w')
    d23sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in d23nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'd23nessaccessorygenesfile' in locals():
            d23nessaccessorygenesfile.close()
        d23nessaccessorygenesfile = open(d23nessaccessorydir + '/' + name + '.fasta', 'w')
    d23nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1esscoregenesfile' in locals():
            st1esscoregenesfile.close()
        st1esscoregenesfile = open(st1esscoredir + '/' + name + '.fasta', 'w')
    st1esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1sesscoregenesfile' in locals():
            st1sesscoregenesfile.close()
        st1sesscoregenesfile = open(st1sesscoredir + '/' + name + '.fasta', 'w')
    st1sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1nesscoregenesfile' in locals():
            st1nesscoregenesfile.close()
        st1nesscoregenesfile = open(st1nesscoredir + '/' + name + '.fasta', 'w')
    st1nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1essaccessorygenesfile' in locals():
            st1essaccessorygenesfile.close()
        st1essaccessorygenesfile = open(st1essaccessorydir + '/' + name + '.fasta', 'w')
    st1essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1sessaccessorygenesfile' in locals():
            st1sessaccessorygenesfile.close()
        st1sessaccessorygenesfile = open(st1sessaccessorydir + '/' + name + '.fasta', 'w')
    st1sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in st1nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'st1nessaccessorygenesfile' in locals():
            st1nessaccessorygenesfile.close()
        st1nessaccessorygenesfile = open(st1nessaccessorydir + '/' + name + '.fasta', 'w')
    st1nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1esscoregenesfile' in locals():
            cs1esscoregenesfile.close()
        cs1esscoregenesfile = open(cs1esscoredir + '/' + name + '.fasta', 'w')
    cs1esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1sesscoregenesfile' in locals():
            cs1sesscoregenesfile.close()
        cs1sesscoregenesfile = open(cs1sesscoredir + '/' + name + '.fasta', 'w')
    cs1sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1nesscoregenesfile' in locals():
            cs1nesscoregenesfile.close()
        cs1nesscoregenesfile = open(cs1nesscoredir + '/' + name + '.fasta', 'w')
    cs1nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1essaccessorygenesfile' in locals():
            cs1essaccessorygenesfile.close()
        cs1essaccessorygenesfile = open(cs1essaccessorydir + '/' + name + '.fasta', 'w')
    cs1essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1sessaccessorygenesfile' in locals():
            cs1sessaccessorygenesfile.close()
        cs1sessaccessorygenesfile = open(cs1sessaccessorydir + '/' + name + '.fasta', 'w')
    cs1sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in cs1nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'cs1nessaccessorygenesfile' in locals():
            cs1nessaccessorygenesfile.close()
        cs1nessaccessorygenesfile = open(cs1nessaccessorydir + '/' + name + '.fasta', 'w')
    cs1nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10esscoregenesfile' in locals():
            h10esscoregenesfile.close()
        h10esscoregenesfile = open(h10esscoredir + '/' + name + '.fasta', 'w')
    salesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10sesscoregenesfile' in locals():
            h10sesscoregenesfile.close()
        h10sesscoregenesfile = open(h10sesscoredir + '/' + name + '.fasta', 'w')
    h10sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10nesscoregenesfile' in locals():
            h10nesscoregenesfile.close()
        h10nesscoregenesfile = open(h10nesscoredir + '/' + name + '.fasta', 'w')
    h10nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10essaccessorygenesfile' in locals():
            h10essaccessorygenesfile.close()
        h10essaccessorygenesfile = open(h10essaccessorydir + '/' + name + '.fasta', 'w')
    h10essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10sessaccessorygenesfile' in locals():
            h10sessaccessorygenesfile.close()
        h10sessaccessorygenesfile = open(h10sessaccessorydir + '/' + name + '.fasta', 'w')
    h10sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in h10nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'h10nessaccessorygenesfile' in locals():
            h10nessaccessorygenesfile.close()
        h10nessaccessorygenesfile = open(h10nessaccessorydir + '/' + name + '.fasta', 'w')
    h10nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2esscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2esscoregenesfile' in locals():
            rh2esscoregenesfile.close()
        rh2esscoregenesfile = open(rh2esscoredir + '/' + name + '.fasta', 'w')
    rh2esscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2sesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2sesscoregenesfile' in locals():
            rh2sesscoregenesfile.close()
        rh2sesscoregenesfile = open(rh2sesscoredir + '/' + name + '.fasta', 'w')
    rh2sesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2nesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2nesscoregenesfile' in locals():
            rh2nesscoregenesfile.close()
        rh2nesscoregenesfile = open(rh2nesscoredir + '/' + name + '.fasta', 'w')
    rh2nesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2essaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2essaccessorygenesfile' in locals():
            rh2essaccessorygenesfile.close()
        rh2essaccessorygenesfile = open(rh2essaccessorydir + '/' + name + '.fasta', 'w')
    rh2essaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2sessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2sessaccessorygenesfile' in locals():
            rh2sessaccessorygenesfile.close()
        rh2sessaccessorygenesfile = open(rh2sessaccessorydir + '/' + name + '.fasta', 'w')
    rh2sessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in rh2nessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'rh2nessaccessorygenesfile' in locals():
            rh2nessaccessorygenesfile.close()
        rh2nessaccessorygenesfile = open(rh2nessaccessorydir + '/' + name + '.fasta', 'w')
    rh2nessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclesscoregenesfile' in locals():
            eclesscoregenesfile.close()
        eclesscoregenesfile = open(eclesscoredir + '/' + name + '.fasta', 'w')
    eclesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclsesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclsesscoregenesfile' in locals():
            eclsesscoregenesfile.close()
        eclsesscoregenesfile = open(eclsesscoredir + '/' + name + '.fasta', 'w')
    eclsesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclnesscoregenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclnesscoregenesfile' in locals():
            eclnesscoregenesfile.close()
        eclnesscoregenesfile = open(eclnesscoredir + '/' + name + '.fasta', 'w')
    eclnesscoregenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclessaccessorygenesfile' in locals():
            eclessaccessorygenesfile.close()
        eclessaccessorygenesfile = open(eclessaccessorydir + '/' + name + '.fasta', 'w')
    eclessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclsessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclsessaccessorygenesfile' in locals():
            eclsessaccessorygenesfile.close()
        eclsessaccessorygenesfile = open(eclsessaccessorydir + '/' + name + '.fasta', 'w')
    eclsessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name

prev_name = ''
for item in eclnessaccessorygenes:
    match_result = match('([a-zA-Z0-9]+?)_\S+', item)
    if not match_result:
        match_result = match('([a-zA-Z]+?)\d+\S*', item)
    name = match_result.group(1)
    if name != prev_name:
        if 'eclnessaccessorygenesfile' in locals():
            eclnessaccessorygenesfile.close()
        eclnessaccessorygenesfile = open(eclnessaccessorydir + '/' + name + '.fasta', 'w')
    eclnessaccessorygenesfile.write(sequences[item].format('fasta'))
    prev_name = name