from os import listdir, path, mkdir
from re import match
from shutil import rmtree

def makedir(dirname):
    # print ('Making directory \'{0}\'..'.format(dirname))
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

class Venn(object):
    def __init__(self, species1, species2, clusters):
        self.sp1 = species1
        self.sp2 = species2
        self.clusters = clusters
        self.sp1presp2abs = 0
        self.sp1presp2pre = 0
        self.sp1abssp2pre = 0
        self.sp1esssp2abs = 0
        self.sp1esssp2nes = 0
        self.sp1esssp2ess = 0
        self.sp1nessp2ess = 0
        self.sp1abssp2ess = 0

    def find_values(self):
        list_of_files = listdir(self.clusters)
        for filename in list_of_files:
            sp1pre = 0
            sp2pre = 0
            sp1ess = 0
            sp2ess = 0
            with open(self.clusters + '/' + filename, 'r') as clusterfile:
                for line in clusterfile:
                    cells = line.split()
                    match_result = match('([a-zA-Z0-9]+?)_\S+', cells[1])
                    if not match_result:
                        match_result = match('([a-zA-Z]+?)\d+\S*', cells[1])
                    if match_result:
                        name = match_result.group(1)
                    if name == self.sp1:
                        sp1pre += 1
                        if float(cells[6]) < 0.2:
                            sp1ess += 1
                    if name == self.sp2:
                        sp2pre += 1
                        if float(cells[6]) < 0.2:
                            sp2ess += 1
            if sp1pre and not sp2pre:
                self.sp1presp2abs += 1
                if sp1ess:
                    self.sp1esssp2abs += 1
            elif sp1pre and sp2pre:
                self.sp1presp2pre += 1
                if sp1ess and sp2ess:
                    self.sp1esssp2ess += 1
                elif sp1ess and not sp2ess:
                    self.sp1esssp2nes += 1
                elif not sp1ess and sp2ess:
                    self.sp1nessp2ess += 1
            elif not sp1pre and sp2pre:
                self.sp1abssp2pre += 1
                if sp2ess:
                    self.sp1abssp2ess += 1


clusts = '/home/fatemeh/EnTrI/results/merge-clust-plot/final_clusters'
output = '/home/fatemeh/EnTrI/results/pairwise-venn-diagrams'
makedir(output)
species_names = ["ROD", "NCTC13441", "Ecoli9000q", "CS17", "ETEC", "t", "SEN", "SL1344", "STM", "STMMW", "ENC", "ERS227112", "BN373"]
name_dict = {"ROD": "Citrobacter_rodentium_ICC168", "CS17": "Escherichia_coli_ETEC_CS17", "ENC": "Enterobacter_cloacae_NCTC_9394",
             "Ecoli9000q": "Escherichia_coli_9000", "ETEC": "Escherichia_coli_ETEC_H10407", "NCTC13441": "Escherichia_coli_UPEC_ST131",
             "ERS227112": "Klebsiella_pneumoniae_RH201207", "BN373": "Klebsiella_pneumoniae_Ecl8", "SEN": "Salmonella_enteritidis_P125109",
             "STM": "Salmonella_typhimurium_A130", "SL1344": "Salmonella_typhimurium_SL1344", "STMMW": "Salmonella_typhimurium_D23580",
             "t": "Salmonella_typhi_Ty2", "b": "Escherichia_coli_K-12"}
venn_values = [[Venn(species_names[i], species_names[j], clusts) for i in range(0, len(species_names))]for j in range(0, len(species_names))]
for i in range(0 , len(species_names)):
    for j in range(i, len(species_names)):
        venn_values[i][j].find_values()

with open(output + '/prepre.txt', 'w') as preprefile:
    preprefile.write('present/present')
    for item in species_names:
        preprefile.write('\t' + name_dict[item])
    preprefile.write('\n')
    for i in range(0,len(species_names)):
        preprefile.write(name_dict[species_names[i]]+'\t')
        for j in range(0, len(species_names)):
            if i <= j:
                preprefile.write(str(venn_values[i][j].sp1presp2pre) + '\t')
            else:
                preprefile.write(str(venn_values[j][i].sp1presp2pre) + '\t')
        preprefile.write('\n')

with open(output + '/preabs.txt', 'w') as preabsfile:
    preabsfile.write('present/absent')
    for item in species_names:
        preabsfile.write('\t' + name_dict[item])
    preabsfile.write('\n')
    for i in range(0,len(species_names)):
        preabsfile.write(name_dict[species_names[i]]+'\t')
        for j in range(0, len(species_names)):
            if i <= j:
                preabsfile.write(str(venn_values[i][j].sp1presp2abs) + '\t')
            else:
                preabsfile.write(str(venn_values[j][i].sp1abssp2pre) + '\t')
        preabsfile.write('\n')

with open(output + '/essess.txt', 'w') as essessfile:
    essessfile.write('essential/essential')
    for item in species_names:
        essessfile.write('\t' + name_dict[item])
    essessfile.write('\n')
    for i in range(0,len(species_names)):
        essessfile.write(name_dict[species_names[i]]+'\t')
        for j in range(0, len(species_names)):
            if i <= j:
                essessfile.write(str(venn_values[i][j].sp1esssp2ess) + '\t')
            else:
                essessfile.write(str(venn_values[j][i].sp1esssp2ess) + '\t')
        essessfile.write('\n')

with open(output + '/essabs.txt', 'w') as essabsfile:
    essabsfile.write('essential/absent')
    for item in species_names:
        essabsfile.write('\t' + name_dict[item])
    essabsfile.write('\n')
    for i in range(0,len(species_names)):
        essabsfile.write(name_dict[species_names[i]]+'\t')
        for j in range(0, len(species_names)):
            if i <= j:
                essabsfile.write(str(venn_values[i][j].sp1esssp2abs) + '\t')
            else:
                essabsfile.write(str(venn_values[j][i].sp1abssp2ess) + '\t')
        essabsfile.write('\n')

with open(output + '/essnes.txt', 'w') as essnesfile:
    essnesfile.write('essential/non-essential')
    for item in species_names:
        essnesfile.write('\t' + name_dict[item])
    essnesfile.write('\n')
    for i in range(0,len(species_names)):
        essnesfile.write(name_dict[species_names[i]]+'\t')
        for j in range(0, len(species_names)):
            if i <= j:
                essnesfile.write(str(venn_values[i][j].sp1esssp2nes) + '\t')
            else:
                essnesfile.write(str(venn_values[j][i].sp1nessp2ess) + '\t')
        essnesfile.write('\n')