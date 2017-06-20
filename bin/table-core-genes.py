from Bio import SeqIO
from re import findall


def read_fasta_sequences(seqdir):
    sequences = {}
    with open(seqdir, 'rU') as fasta_file:
        sequences.update(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')))
    return sequences

core_essential_pth = '../results/define-core-accessory-hieranoid/all/core-essential-genomes/b.fasta'
core_sometimes_essential_pth = '../results/define-core-accessory-hieranoid/all/core-sometimes-essential-genomes/b.fasta'
core_never_essential_pth = '../results/define-core-accessory-hieranoid/all/core-never-essential-genomes/b.fasta'

ancestrally_essential_pth = '../results/define-core-accessory-hieranoid-fitch/all/core-essential-genomes/b.fasta'
ancestrally_sometimes_essential_pth = '../results/define-core-accessory-hieranoid-fitch/all/core-sometimes-essential-genomes/b.fasta'
ancestrally_never_essential_pth = '../results/define-core-accessory-hieranoid-fitch/all/core-never-essential-genomes/b.fasta'

seqdb_pth = '../data/fasta-protein/chromosome/U00096.fasta'

pathways_pth = '../results/KEGG/escherichia_coli_K-12_MG1655.dat'
kegg_groups_pth = '../results/table-core-genes/kegg-groups.txt'

outdir = '../results/table-core-genes/table.tsv'

core_essential = read_fasta_sequences(core_essential_pth)
core_sometimes_essential = read_fasta_sequences(core_sometimes_essential_pth)
core_never_essential = read_fasta_sequences(core_never_essential_pth)

ancestrally_essential = read_fasta_sequences(ancestrally_essential_pth)
ancestrally_sometimes_essential = read_fasta_sequences(ancestrally_sometimes_essential_pth)
ancestrally_never_essential = read_fasta_sequences(ancestrally_never_essential_pth)

seqdb = read_fasta_sequences(seqdb_pth)

with open(kegg_groups_pth, 'r') as groupsfile:
    with open(outdir, 'w') as tofile:
        tofile.write('Biological process\tSubprocess\tCore essential\tCore ambiguous\tCore non-essential\tNot core\t' +
                     'Ancestrally essential\tAncestrally ambiguous\tAncestrally non-essential\tNot ancestral\n')
        for gline in groupsfile:
            if not gline.startswith('>'):
                pathway_name = gline[:-1]
                coes = list()
                coso = list()
                cone = list()
                noco = list()
                anes = list()
                anso = list()
                anne = list()
                noan = list()
                with open(pathways_pth, 'r') as pathwaysfile:
                    for pline in pathwaysfile:
                        cells = pline.split('\t')
                        if cells[0] == '\"gene_id\"':
                            continue
                        pathway_name2 = cells[2].split("-")[0]
                        pathway_name2 = pathway_name2[1:-1]
                        if pathway_name == pathway_name2:
                            gene_tag = cells[0][1:-1]
                            if gene_tag not in seqdb:
                                continue
                            desc = seqdb[gene_tag].description
                            gene_name = findall('\[([a-zA-Z]{0,5})\]', desc)[
                                len(findall('\[([a-zA-Z]{0,5})\]', desc)) - 1]
                            if gene_tag in core_essential:
                                coes.append(gene_name)
                            elif gene_tag in core_sometimes_essential:
                                coso.append(gene_name)
                            elif gene_tag in core_never_essential:
                                cone.append(gene_name)
                            else:
                                noco.append(gene_name)

                            if gene_tag in ancestrally_essential:
                                anes.append(gene_name)
                            elif gene_tag in ancestrally_sometimes_essential:
                                anso.append(gene_name)
                            elif gene_tag in ancestrally_never_essential:
                                anne.append(gene_name)
                            else:
                                noan.append(gene_name)

                if coes or coso or cone or noco or anes or anso or anne or noan:
                    coes.sort()
                    coso.sort()
                    cone.sort()
                    noco.sort()
                    anes.sort()
                    anso.sort()
                    anne.sort()
                    noan.sort()

                    i = 0
                    coes2 = list()
                    while i < len(coes)-1:
                        item = coes[i]
                        if len(item) != 4:
                            coes2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(coes) and len(coes[j]) == 4 and coes[j][0:3] == item[0:3]:
                                item2 = coes[j]
                                string += item2[3]
                                j += 1
                            i = j
                            coes2.append(string)

                    i = 0
                    coso2 = list()
                    while i < len(coso) - 1:
                        item = coso[i]
                        if len(item) != 4:
                            coso2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(coso) and len(coso[j]) == 4 and coso[j][0:3] == item[0:3]:
                                item2 = coso[j]
                                string += item2[3]
                                j += 1
                            i = j
                            coso2.append(string)

                    i = 0
                    cone2 = list()
                    while i < len(cone) - 1:
                        item = cone[i]
                        if len(item) != 4:
                            cone2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(cone) and len(cone[j]) == 4 and cone[j][0:3] == item[0:3]:
                                item2 = cone[j]
                                string += item2[3]
                                j += 1
                            i = j
                            cone2.append(string)
                    i = 0
                    noco2 = list()
                    while i < len(noco) - 1:
                        item = noco[i]
                        if len(item) != 4:
                            noco2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(noco) and len(noco[j]) == 4 and noco[j][0:3] == item[0:3]:
                                item2 = noco[j]
                                string += item2[3]
                                j += 1
                            i = j
                            noco2.append(string)

                    i = 0
                    anes2 = list()
                    while i < len(anes) - 1:
                        item = anes[i]
                        if len(item) != 4:
                            anes2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(anes) and len(anes[j]) == 4 and anes[j][0:3] == item[0:3]:
                                item2 = anes[j]
                                string += item2[3]
                                j += 1
                            i = j
                            anes2.append(string)

                    i = 0
                    anso2 = list()
                    while i < len(anso) - 1:
                        item = anso[i]
                        if len(item) != 4:
                            anso2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(anso) and len(anso[j]) == 4 and anso[j][0:3] == item[0:3]:
                                item2 = anso[j]
                                string += item2[3]
                                j += 1
                            i = j
                            anso2.append(string)

                    i = 0
                    anne2 = list()
                    while i < len(anne) - 1:
                        item = anne[i]
                        if len(item) != 4:
                            anne2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(anne) and len(anne[j]) == 4 and anne[j][0:3] == item[0:3]:
                                item2 = anne[j]
                                string += item2[3]
                                j += 1
                            i = j
                            anne2.append(string)

                    i = 0
                    noan2 = list()
                    while i < len(noan) - 1:
                        item = noan[i]
                        if len(item) != 4:
                            noan2.append(item)
                            i += 1
                        else:
                            string = item
                            j = i + 1
                            while j < len(noan) and len(noan[j]) == 4 and noan[j][0:3] == item[0:3]:
                                item2 = noan[j]
                                string += item2[3]
                                j += 1
                            i = j
                            noan2.append(string)

                    coes_print = str(coes2)[1:-1].replace('\'', '')
                    coso_print = str(coso2)[1:-1].replace('\'', '')
                    cone_print = str(cone2)[1:-1].replace('\'', '')
                    noco_print = str(noco2)[1:-1].replace('\'', '')
                    anes_print = str(anes2)[1:-1].replace('\'', '')
                    anso_print = str(anso2)[1:-1].replace('\'', '')
                    anne_print = str(anne2)[1:-1].replace('\'', '')
                    noan_print = str(noan2)[1:-1].replace('\'', '')
                    tofile.write(process + '\t')
                    tofile.write(pathway_name + '\t')
                    tofile.write(coes_print + '\t')
                    tofile.write(coso_print + '\t')
                    tofile.write(cone_print + '\t')
                    tofile.write(noco_print + '\t')
                    tofile.write(anes_print + '\t')
                    tofile.write(anso_print + '\t')
                    tofile.write(anne_print + '\t')
                    tofile.write(noan_print + '\n')
            else:
                process = gline[1:-1]