from collections import defaultdict
from re import match
from os import listdir, path
from math import floor

seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
plots = '../data/plot-files/chromosome'
essentialitydir = '../results/biases/dbscan/'
outdir = '../results/5primeinserts.txt'
lfiles = listdir(plots)
lfiles.remove('BW25113.plot')
lfiles.remove('EC958.plot')
list_of_list_of_files = [lfiles, ['BW25113.plot'], ['EC958.plot']]
# list_of_list_of_files = [[i] for i in lfiles]
for item in list_of_list_of_files:
    forward_dict = defaultdict(list)
    reverse_dict = defaultdict(list)
    for filename in item:
        name, extension = path.splitext(filename)
        with open('{0}/{1}'.format(plots, filename), 'r') as plotfile:
            for line in plotfile:
                cells = line.split()
                forward_dict[name].append(int(int(cells[0]) > 0))
                reverse_dict[name].append(int(int(cells[1]) > 0))

    list_of_files = listdir(essentialitydir)
    essentiality_dict = dict()
    for filename in list_of_files:
        with open(essentialitydir + filename) as fromfile:
            for line in fromfile:
                cells = line.split()
                essentiality_dict[cells[0]] = cells[2]

    ff = 0  # gene: forward, transposon: forward
    ff0 = 0  # frame: +0
    ff1 = 0
    ff2 = 0
    fr = 0  # gene: forward, transposon: reverse
    fr0 = 0
    fr1 = 0
    fr2 = 0
    rf = 0  # gene: reverse, transposon: forward
    rf0 = 0
    rf1 = 0
    rf2 = 0
    rr = 0  # gene: reverse, transposon: reverse
    rr0 = 0
    rr1 = 0
    rr2 = 0
    counter = 0
    with open(seqdb, 'r') as sequencefile:
        for line in sequencefile:
            if line.startswith('>'):
                match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+)\s+\[\S+/(\d+)\-(\d+)\s\((\w+)\)', line)
                if match_result is not None:
                    gene_name = match_result.group(1)
                    strain_name = match_result.group(2)
                    start = int(match_result.group(3)) - 1
                    end = int(match_result.group(4)) - 1
                    strand = match_result.group(5)

                    if strain_name in forward_dict.keys() and gene_name in essentiality_dict.keys() \
                            and essentiality_dict[gene_name] == 'essential':
                        if strand == 'Forward':
                            fiveperc = round((end - start + 1) * 5 / 100) + start
                            if sum(forward_dict[strain_name][start:fiveperc + 1]) > 0 or sum(
                                    reverse_dict[strain_name][start:fiveperc + 1]) > 0:
                                counter += 1
                            ff += sum(forward_dict[strain_name][start:fiveperc + 1])
                            ff0 += sum(
                                [forward_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 0])
                            ff1 += sum(
                                [forward_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 1])
                            ff2 += sum(
                                [forward_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 2])
                            fr += sum(reverse_dict[strain_name][start:fiveperc + 1])
                            fr0 += sum(
                                [reverse_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 0])
                            fr1 += sum(
                                [reverse_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 1])
                            fr2 += sum(
                                [reverse_dict[strain_name][i] for i in range(start, fiveperc + 1) if
                                 (i - start) % 3 == 2])
                        elif strand == 'Complement':
                            fiveperc = end - round((end - start + 1) * 5 / 100)
                            if sum(forward_dict[strain_name][fiveperc:end + 1]) > 0 or sum(
                                    reverse_dict[strain_name][fiveperc:end + 1]) > 0:
                                counter += 1
                            rf += sum(forward_dict[strain_name][fiveperc:end + 1])
                            rf0 += sum(
                                [forward_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 0])
                            rf1 += sum(
                                [forward_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 1])
                            rf2 += sum(
                                [forward_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 2])
                            rr += sum(reverse_dict[strain_name][fiveperc:end + 1])
                            rr0 += sum(
                                [reverse_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 0])
                            rr1 += sum(
                                [reverse_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 1])
                            rr2 += sum(
                                [reverse_dict[strain_name][i] for i in range(fiveperc, end + 1) if
                                 (end - i) % 3 == 2])

    print(item)
    print('1: ' + str(ff0 + rr0) + ' = ' + str(ff0) + ' + ' + str(rr0))
    print('2: ' + str(ff1 + rr1) + ' = ' + str(ff1) + ' + ' + str(rr1))
    print('3: ' + str(ff2 + rr2) + ' = ' + str(ff2) + ' + ' + str(rr2))
    print('4: ' + str(fr0 + rf0) + ' = ' + str(fr0) + ' + ' + str(rf0))
    print('5: ' + str(fr1 + rf1) + ' = ' + str(fr1) + ' + ' + str(rf1))
    print('6: ' + str(fr2 + rf2) + ' = ' + str(fr2) + ' + ' + str(rf2))
    print('\n')

    if item == lfiles:
        print(str(counter)+'\n')
        with open(outdir, 'w') as tofile:
            tofile.write(str(ff0 + rr0) + '\n' + str(ff1 + rr1) + '\n' + str(ff2 + rr2) + '\n')
            tofile.write(str(fr0 + rf0) + '\n' + str(fr1 + rf1) + '\n' + str(fr2 + rf2) + '\n')
    # print('forward gene, forward transposon: ' + str(ff))
    # print('+0: ' + str(ff0))
    # print('+1: ' + str(ff1))
    # print('+2: ' + str(ff2))
    # print('forward gene, reverse transposon: ' + str(fr))
    # print('+0: ' + str(fr0))
    # print('+1: ' + str(fr1))
    # print('+2: ' + str(fr2))
    # print('reverse gene, forward transposon: ' + str(rf))
    # print('+0: ' + str(rf0))
    # print('+1: ' + str(rf1))
    # print('+2: ' + str(rf2))
    # print('reverse gene, reverse transposon: ' + str(rr))
    # print('+0: ' + str(rr0))
    # print('+1: ' + str(rr1))
    # print('+2: ' + str(rr2))