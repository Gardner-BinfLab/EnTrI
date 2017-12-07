from os import listdir

def read_essentiality(filespath):
    list_of_files = listdir(filespath)
    ess_dict = {}
    for filename in list_of_files:
        with open(filespath + filename) as fromfile:
            for line in fromfile:
                cells = line.split()
                if cells[2] in ['essential', 'non-essential', 'beneficial-loss']:
                    ess_dict[cells[0]] = cells[2]
    return ess_dict

essentiality = '../results/biases/dbscan/'
ess_dict = read_essentiality(essentiality)
hmmscan = '../results/dcgo/seqdb_Pfam-A.domtblout'
essen = '../results/dcgo/essential_pfam-A.domtblout'
noness = '../results/dcgo/non-essential_pfam-A.domtblout'
benloss = '../results/dcgo/beneficial-loss_pfam-A.domtblout'

with open(hmmscan, 'r') as fromfile:
    with open(essen, 'w') as essenfile:
        with open(noness, 'w') as nonessfile:
            with open(benloss, 'w') as benlossfile:
                for line in fromfile:
                    if not line.startswith('#'):
                        cells = line.split()
                        if cells[3] in ess_dict:
                            if ess_dict[cells[3]] == 'essential':
                                essenfile.write(line)
                            elif ess_dict[cells[3]] == 'non-essential':
                                nonessfile.write(line)
                            else:
                                benlossfile.write(line)


