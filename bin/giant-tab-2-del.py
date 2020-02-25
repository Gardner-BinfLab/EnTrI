import pandas as pd

outdir = '../results/giant-tab/'
enteropath = outdir + 'entero_annotated_clusters.txt'
degpath = outdir + 'deg_annotated_clusters.txt'
endopath = outdir + 'endo_annotated_clusters.txt'
outpath = outdir + 'giant-table.tsv'

enteroannot = pd.read_csv(enteropath, sep='\t', header=None, names=['bactNOG','gene','c1','c2','c3','c4','c5','c6','c7',
                                                                    'c8','c9','c10','c11','c12','c13'])
enteroannot.set_index(['bactNOG','gene'],inplace=True)
degannot = pd.read_csv(enteropath, sep='\t', header=None, names=['bactNOG','gene','c1','c2','c3','c4','c5','c6','c7',
                                                                 'c8','c9','c10','c11','c12','c13','c14','c15','c16',
                                                                 'c17','c18','c19','c20','c21','c22','c23','c24','c25',
                                                                 'c26','c27'])
degannot.set_index(['bactNOG','gene'],inplace=True)
endoannot = pd.read_csv(enteropath, sep='\t', header=None, names=['bactNOG','gene','c1','c2','c3','c4','c5','c6','c7',
                                                                 'c8','c9','c10','c11','c12','c13','c14','c15','c16',
                                                                 'c17','c18','c19','c20','c21','c22','c23','c24','c25',
                                                                 'c26','c27','c28','c29','c30','c31','c32','c33','c34'])
endoannot.set_index(['bactNOG','gene'],inplace=True)

colnames = [
    'Citrobacter rodentium ICC168', 'Salmonella Typhi Ty2', 'Salmonella Enteritidis P125109',
    'Salmonella Typhimurium SL1344', 'Salmonella Typhimurium SL3261', 'Salmonella Typhimurium D23580',
    'Salmonella Typhimurium A130', 'Escherichia coli UPEC ST131 NCTC13441', 'Escherichia coli ST131 EC958',
    'Escherichia coli BW25113', 'Escherichia coli BW25113 (Keio)', 'Klebsiella pneumoniae RH201207',
    'Klebsiella pneumoniae Ecl8',
    'Mycobacterium tuberculosis H37Rv', 'Synechococcus elongatus PCC 7942', 'Porphyromonas gingivalis ATCC 33277',
    'Bacteroides thetaiotaomicron VPI-5482','Bacteroides fragilis 638R', 'Caulobacter crescentus',
    'Brevundimonas subvibrioides ATCC 15264', 'Rhodopseudomonas palustris CGA009', 'Agrobacterium fabrum str. C58',
    'Burkholderia thailandensis E264', 'Burkholderia pseudomallei K96243', 'Francisella novicida U112',
    'Shewanella oneidensis MR-1', 'Vibrio cholerae N16961', 'Haemophilus influenzae Rd KW20',
    'Salmonella enterica serovar Typhi', 'Acinetobacter baumannii ATCC 17978',
    'Pseudomonas aeruginosa PAO1', 'Helicobacter pylori 26695', 'Mycoplasma genitalium G37',
    'Mycoplasma pulmonis UAB CTIP', 'Streptococcus agalactiae A909', 'Streptococcus pyogenes MGAS5448',
    'Streptococcus pyogenes NZ131', 'Staphylococcus aureus NCTC 8325', 'Staphylococcus aureus N315',
    'Bacillus subtilis 168', 'Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
    'Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans', 'Sodalis glossinidius str. morsitans',
    'Sodalis praecaptivus strain HS1', 'Candidatus Sodalis pierantonius str. SOPE', 'Secondary endosymbiont of Ctenarytaina eucalypti',
    'Candidatus Moranella endobia PCIT', 'Secondary endosymbiont of Heteropsylla cubana',
    'Candidatus Baumannia cicadellinicola strain BGSS', 'Candidatus Baumannia cicadellinicola strain B-GSS',
    'Baumannia cicadellinicola str. Hc (Homalodisca coagulata)', 'Candidatus Blochmannia chromaiodes str. 640',
    'Candidatus Blochmannia pennsylvanicus str. BPEN', 'Candidatus Blochmannia floridanus', 'Candidatus Blochmannia vafer str. BVAF',
    'Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675',
    'Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757', 'Buchnera aphidicola str. Bp (Baizongia pistaciae)',
    'Buchnera aphidicola BCc', 'Buchnera aphidicola (Cinara tujafilina)', 'Buchnera aphidicola str. F009 (Myzus persicae)',
    'Buchnera aphidicola str. W106 (Myzus persicae)', 'Buchnera aphidicola str. USDA (Myzus persicae)',
    'Buchnera aphidicola str. G002 (Myzus persicae)', 'Buchnera aphidicola str. Ua (Uroleucon ambrosiae)',
    'Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)', 'Buchnera aphidicola str. APS (Acyrthosiphon pisum)',
    'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)',
    'Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)',
    'Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. 5A (Acyrthosiphon pisum)',
    'Buchnera aphidicola str. Sg (Schizaphis graminum)']
tags = ['ROD', 't', 'SEN', 'SL1344', 'SL3261', 'STMMW', 'STM', 'NCTC13441', 'EC958', 'BW25113', 'b', 'ERS227112',
        'BN373', 'DEG1027', 'DEG1040', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1020', 'DEG1046', 'DEG1041', 'DEG1045',
        'DEG1024', 'DEG1035', 'DEG1012', 'DEG1029', 'DEG1003', 'DEG1005', 'DEG1016', 'DEG1043', 'DEG1036', 'DEG1008',
        'DEG1006', 'DEG1014', 'DEG1042', 'DEG1037', 'DEG1038', 'DEG1017', 'DEG1002', 'DEG1001', 'BA000021', 'WIGMOR',
        'SG', 'Sant', 'SOPEG', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BCHRO640', 'BPEN', 'Bfl', 'BVAF',
        'BTURN675', 'BOBLI757', 'bbp', 'BCc', 'BCTU', 'BUMPF009', 'BUMPW106', 'BUMPUSDA', 'BUMPG002', 'BUAMB', 'BAKON',
        'BA000003', 'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUsg']

with open(outpath, 'w') as tofile:
    for item in colnames[:-1]:
        tofile.write(item + '\t')
    tofile.write(colnames[-1] + '\n')

    for index, row in enteroannot.iterrows():
        dict = {el: 0 for el in tags}




