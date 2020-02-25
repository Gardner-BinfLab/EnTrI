from os import listdir, path, mkdir
from shutil import rmtree
from re import match
import pandas as pd

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

def read_k12(inpath):
    all = []
    with open(inpath) as from_file:
        for line in from_file:
            cell = line.split()
            all.append(cell[0])
    return all

outdir = '../results/giant-tab/'
makedir(outdir)
outpath = outdir + 'giant-tab.tsv'

colnames = ['DEG: Synechococcus elongatus PCC 7942', 'DEG: Mycobacterium tuberculosis H37Rv',
            'DEG: Porphyromonas gingivalis ATCC 33277', 'DEG: Bacteroides thetaiotaomicron VPI-5482',
            'DEG: Bacteroides fragilis 638R', 'DEG: Francisella novicida U112', 'DEG: Acinetobacter baumannii ATCC 17978',
            'DEG: Pseudomonas aeruginosa PAO1', 'DEG: Shewanella oneidensis MR-1', 'TraDIS: Klebsiella pneumoniae Ecl8',
            'TraDIS: Klebsiella pneumoniae RH201207', 'TraDIS: Escherichia coli ST131 EC958', 'TraDIS: Escherichia coli UPEC ST131 NCTC13441',
            'Keio: Escherichia coli BW25113', 'TraDIS: Escherichia coli BW25113', 'TraDIS: Citrobacter rodentium ICC168',
            'TraDIS: Salmonella Typhi Ty2', 'DEG: Salmonella enterica serovar Typhi', 'TraDIS: Salmonella Typhimurium A130',
            'TraDIS: Salmonella Typhimurium D23580', 'TraDIS: Salmonella Typhimurium SL3261', 'TraDIS: Salmonella Typhimurium SL1344',
            'TraDIS: Salmonella Enteritidis P125109', 'Symbiont: Secondary endosymbiont of Ctenarytaina eucalypti',
            'Symbiont: Candidatus Moranella endobia PCIT', 'Symbiont: Secondary endosymbiont of Heteropsylla cubana',
            'Symbiont: Candidatus Baumannia cicadellinicola strain BGSS', 'Symbiont: Candidatus Baumannia cicadellinicola strain B-GSS',
            'Symbiont: Baumannia cicadellinicola str. Hc (Homalodisca coagulata)',
            'Symbiont: Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757',
            'Symbiont: Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675',
            'Symbiont: Candidatus Blochmannia vafer str. BVAF', 'Symbiont: Candidatus Blochmannia floridanus',
            'Symbiont: Candidatus Blochmannia pennsylvanicus str. BPEN', 'Symbiont: Candidatus Blochmannia chromaiodes str. 640',
            'Symbiont: Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
            'Symbiont: Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans',
            'Symbiont: Buchnera aphidicola str. Sg (Schizaphis graminum)', 'Symbiont: Buchnera aphidicola str. G002 (Myzus persicae)',
            'Symbiont: Buchnera aphidicola str. USDA (Myzus persicae)', 'Symbiont: Buchnera aphidicola str. F009 (Myzus persicae)',
            'Symbiont: Buchnera aphidicola str. W106 (Myzus persicae)', 'Symbiont: Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)',
            'Symbiont: Buchnera aphidicola str. APS (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)', 'Symbiont: Buchnera aphidicola str. 5A (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. Ua (Uroleucon ambrosiae)', 'Symbiont: Buchnera aphidicola str. Bp (Baizongia pistaciae)',
            'Symbiont: Buchnera aphidicola BCc', 'Symbiont: Buchnera aphidicola (Cinara tujafilina)',
            'Symbiont: Sodalis glossinidius str. morsitans', 'Symbiont: Sodalis praecaptivus strain HS1',
            'Symbiont: Candidatus Sodalis pierantonius str. SOPE', 'DEG: Haemophilus influenzae Rd KW20', 'DEG: Vibrio cholerae N16961',
            'DEG: Burkholderia pseudomallei K96243', 'DEG: Burkholderia thailandensis E264', 'DEG: Rhodopseudomonas palustris CGA009',
            'DEG: Agrobacterium fabrum str. C58', 'DEG: Caulobacter crescentus', 'DEG: Brevundimonas subvibrioides ATCC 15264',
            'DEG: Helicobacter pylori 26695', 'DEG: Mycoplasma genitalium G37', 'DEG: Mycoplasma pulmonis UAB CTIP',
            'DEG: Streptococcus agalactiae A909', 'DEG: Streptococcus pyogenes NZ131', 'DEG: Streptococcus pyogenes MGAS5448',
            'DEG: Staphylococcus aureus NCTC 8325', 'DEG: Staphylococcus aureus N315', 'DEG: Bacillus subtilis 168',
            'Locus: Klebsiella pneumoniae Ecl8',
            'Locus: Klebsiella pneumoniae RH201207', 'Locus: Escherichia coli ST131 EC958',
            'Locus: Escherichia coli UPEC ST131 NCTC13441',
            'Locus: Escherichia coli BW25113 (Keio)', 'Locus: Escherichia coli BW25113',
            'Locus: Citrobacter rodentium ICC168',
            'Locus: Salmonella Typhi Ty2',
            'Locus: Salmonella Typhimurium A130',
            'Locus: Salmonella Typhimurium D23580', 'Locus: Salmonella Typhimurium SL3261',
            'Locus: Salmonella Typhimurium SL1344',
            'Locus: Salmonella Enteritidis P125109'
            ]

tags = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008', 'DEG1006',
        'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002', 'DEG1001']
genome_dict = dict(zip(tags, colnames))
tradis = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'STM', 'STMMW', 'SL3261',
                'SL1344', 'SEN']
enterobacter = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
                'SL3261', 'SL1344', 'SEN']
endosymbiont = ['A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
                'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON',
                'BA000003', 'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG',
                'Sant', 'SOPEG']
gammaprot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003']
prot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
gammaprotexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b',
                   'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003']
protexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113',
              'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035',
              'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
bactexsymb = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029',
              'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
              'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020',
              'DEG1046', 'DEG1008', 'DEG1006', 'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002',
              'DEG1001']
lenenterobacter = len(enterobacter)
lenendosymbiont = len(endosymbiont)
lengammaprot = len(gammaprot)
lenprot = len(prot)
lengammaprotexsymb = len(gammaprotexsymb)
lenprotexsymb = len(protexsymb)
lenbact = len(tags)
lenbactexsymb = len(bactexsymb)

essentiality = '../results/biases/dbscan/'
ii_dict = dict()
list_of_files = listdir(essentiality)
for filename in list_of_files:
    with open(essentiality + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            ii_dict[cells[0]] = cells[3]

annotations = '../results/all-hieranoid/cluster-representatives.emapper.annotations'
annot_dict = {}
with open(annotations, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('#'):
            cells = line.split('\t')
            annot_dict[cells[0]] = cells[5]
list_annot_keys = list(annot_dict.keys())

k12path = '../results/ecogene-k12.txt'
k12genes = read_k12(k12path)

kegg = '../results/KEGG/escherichia_coli_K-12_MG1655.dat'
kegg_dict = {}
with open(kegg, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[2] = cells[2][1:-2]
        if cells[0] != 'gene_id' and cells[0] not in kegg_dict.keys():
            cells[2] = match('(.*) - Escherichia coli K-12 MG1655', cells[2]).group(1)
            kegg_dict[cells[0]] = cells[2]
        elif cells[0] != 'gene_id':
            cells[2] = match('(.*) - Escherichia coli K-12 MG1655', cells[2]).group(1)
            kegg_dict[cells[0]] += ' / ' + cells[2]

clusters = '../results/all-hieranoid/clusters.txt'
prev_name = ''
with open(outpath, 'w') as tofile:
    tofile.write('Gene\tPathway\t')
    for item in colnames:
        tofile.write(item + '\t')
    tofile.write('Enterobacteriaceae %essential\tEndosymbiont %conserved\tGammaproteobacteria %essential\t' +
                 'Gammaproteobacteria (excluding symbionts) %essential\tProteobacteria %essential\t' +
                 'Proteobacteria (excluding symbionts) %essential\tBacteria %essential\t' +
                 'Bacteria (excluding symbionts) %essential\n'
                 )
    with open(clusters, 'r') as fromfile:
        for line in fromfile:
            pathway = '-'
            genes = line.split()
            genes = [item for item in genes if not item.startswith('exDEG')]
            intersect = list(set(genes) & set(list_annot_keys))
            if len(intersect) == 1:
                intersect = intersect[0]
                if len(genes) > 2:
                    if annot_dict[intersect]:
                        clustdict = {}
                        for item in tags:
                            if item.startswith('DEG'):
                                clustdict[item] = 'A/N'
                            else:
                                clustdict[item] = 'A'
                        locusdict = {el:'-' for el in tradis}

                        numenterobacter = 0
                        numendosymbiont = 0
                        numgammaprot = 0
                        numprot = 0
                        numgammaprotexsymb = 0
                        numprotexsymb = 0
                        numbact = 0
                        numbactexsymb = 0

                        tofile.write(annot_dict[intersect] + '\t')
                        for gene in genes:
                            if gene.startswith('DEG'):
                                species = gene[0:7]
                            else:
                                species = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', str(gene)).group(1).strip('_')
                            if species in tradis:
                                locusdict[species] = gene
                                if species != 'b':
                                    if gene in ii_dict.keys():
                                        clustdict[species] = str(ii_dict[gene])
                                        if float(ii_dict[gene]) <= 0:
                                            if species in enterobacter:
                                                numenterobacter += 1
                                            if species in endosymbiont:
                                                numendosymbiont += 1
                                            if species in gammaprot:
                                                numgammaprot += 1
                                            if species in gammaprotexsymb:
                                                numgammaprotexsymb += 1
                                            if species in prot:
                                                numprot += 1
                                            if species in protexsymb:
                                                numprotexsymb += 1
                                            if species in tags:
                                                numbact += 1
                                            if species in bactexsymb:
                                                numbactexsymb += 1
                                else:
                                    if gene in k12genes:
                                        clustdict[species] = 'E'
                                        if species in enterobacter:
                                            numenterobacter += 1
                                        if species in endosymbiont:
                                            numendosymbiont += 1
                                        if species in gammaprot:
                                            numgammaprot += 1
                                        if species in gammaprotexsymb:
                                            numgammaprotexsymb += 1
                                        if species in prot:
                                            numprot += 1
                                        if species in protexsymb:
                                            numprotexsymb += 1
                                        if species in tags:
                                            numbact += 1
                                        if species in bactexsymb:
                                            numbactexsymb += 1
                                    else:
                                        clustdict[species] = 'N'
                                    if gene in kegg_dict.keys():
                                        pathway = kegg_dict[gene]
                            else:
                                if gene.startswith('DEG'):
                                    clustdict[species] = 'E'
                                else:
                                    clustdict[species] = 'P'

                                if species in enterobacter:
                                    numenterobacter += 1
                                if species in endosymbiont:
                                    numendosymbiont += 1
                                if species in gammaprot:
                                    numgammaprot += 1
                                if species in gammaprotexsymb:
                                    numgammaprotexsymb += 1
                                if species in prot:
                                    numprot += 1
                                if species in protexsymb:
                                    numprotexsymb += 1
                                if species in tags:
                                    numbact += 1
                                if species in bactexsymb:
                                    numbactexsymb += 1

                        tofile.write(pathway + '\t')
                        for item in tags:
                            tofile.write(clustdict[item] + '\t')
                        for item in tradis:
                            tofile.write(locusdict[item] + '\t')
                        percenterobacter = round(numenterobacter*100 / lenenterobacter,1)
                        percendosymbiont = round(numendosymbiont * 100 / lenendosymbiont, 1)
                        percgammaprot = round(numgammaprot * 100 / lengammaprot, 1)
                        percgammaprotexsymb = round(numgammaprotexsymb * 100 / lengammaprotexsymb, 1)
                        percprot = round(numprot * 100 / lenprot, 1)
                        percprotexsymb = round(numprotexsymb * 100 / lenprotexsymb, 1)
                        percbact = round(numbact * 100 / lenbact, 1)
                        percbactexsymb = round(numbactexsymb * 100 / lenbactexsymb, 1)
                        tofile.write(str(percenterobacter) + '\t' + str(percendosymbiont) + '\t' + str(percgammaprot) +
                                     '\t' + str(percgammaprotexsymb) + '\t' + str(percprot) + '\t' +
                                     str(percprotexsymb) + '\t' + str(percbact) + '\t' + str(percbactexsymb) + '\n')

                        prev_name = annot_dict[intersect]

tbl = pd.read_csv(outpath, sep='\t', header=0)
tbl = tbl.iloc[tbl.Gene.str.lower().argsort()]
tbl.reset_index(inplace=True)
tbl.drop('index', axis=1, inplace=True)
tbl.to_csv(outpath, sep='\t', index=False)