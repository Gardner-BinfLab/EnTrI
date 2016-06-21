from collections import defaultdict
from re import match, findall
from os import listdir, path

seqdb = '../sequences/fasta-dna/chromosome/seqdb.fasta'
iis = '../results/insertion-indices/gamma'
results = '../results/insertion-indices/check-biases'
genome_length = {"SL1344":4878012, "STMMW":4879400, "SEN":4685848, "t":4791961, "STM":4895639, "ETEC":5153435,
                 "b":4641652, "CS17":4994793, "NCTC13441":5174631, "ROD":5346659, "BN373":5324709, "ERS227112":5869288,
                 "ENC":4908759}

list_of_files = listdir(iis)
for filename in list_of_files:
    iis_dict = defaultdict(tuple)
    with open('{0}/{1}'.format(iis, filename), 'r') as iifile:
        for line in iifile:
            cells = line.split()
            iis_dict[cells[0]] = (float(cells[1]), cells[2])

    gene_name = ''
    with open(results + '/' + filename, 'w') as tofile:
        with open(seqdb, 'r') as sequencefile:
            for line in sequencefile:
                if line.startswith('>'):
                    if gene_name != '':
                        tofile.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(gene_name, iis_dict[gene_name][0], iis_dict[gene_name][1], float(start[0])/genome_length[strain_name], float(gc)/gene_length))
                    gc = 0
                    match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/((\d+\-\d+\s)+)\(', line)
                    if match_result is None:
                        match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/((\d+\-\d+\s)+)\(', line)
                    if match_result is not None:
                        gene_name = match_result.group(1)
                        strain_name = match_result.group(2)
                        starts_ends = findall('\d+\-\d+\s', match_result.group(3))
                        start = []
                        end =  []
                        for item in starts_ends:
                            match_result = match('(\d+)-(\d+)', item)
                            start.append(int(match_result.group(1)))
                            end.append(int(match_result.group(2)))
                    else:
                        gene_name = ''
                    if gene_name not in iis_dict.keys():
                        gene_name = ''
                    else:
                        # gene_insertions = 0
                        gene_length = 0
                        for i in range(0,len(start)):
                            # gene_insertions += sum([1 for x in iis_dict[strain_name][start[i]-1:end[i]] if x > 0])
                            gene_length += end[i] - start[i] + 1
                        # ii = (float(gene_insertions)/gene_length)/(float(genome_insertions[strain_name])/genome_length[strain_name])
                else:
                    gc += line.count('g') + line.count('c')