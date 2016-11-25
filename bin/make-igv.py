from re import match

fastafile = '../data/fasta-protein/chromosome/CS17.fasta'
plotfile = '../data/plot-files/chromosome/CS17.plot'
igvfile = '../results/igvs/CS17.igv'

plots = []
with open(plotfile, 'r') as fromfile:
    for line in fromfile:
        cells = line.split()
        plots.append(int(cells[0])+int(cells[1]))

with open(fastafile, 'r') as fromfile:
    with open(igvfile, 'w') as tofile:
        tofile.write('#type=COPY_NUMBER\nChromosome\tStart\tEnd\tFeature\tReads\tTAs\n')
        tas = str(1)
        for line in fromfile:
            if line.startswith('>'):
                #match_result = match('>(\S+)\s+\[[^/]+/(\d+)\-(\d+)', line)
                match_result = match('>\s*((\S+?)_+\S+)\s+\[[^/]+/(\d+)\-(\d+)', line)
                if match_result is None:
                    match_result = match('>\s*(([a-zA-Z]+)\d+[a-zA-Z]?)\s+\[[^/]+/(\d+)\-(\d+)', line)
                if match_result is not None:
                    feature = match_result.group(1)
                    chromosome = match_result.group(2)
                    genestart = int(match_result.group(3))
                    geneend = int(match_result.group(4))
                    for i in range(genestart-1, geneend):
                        reads = str(plots[i])
                        start = str(i)
                        end = str(i)
                        tofile.write(chromosome + '\t' + start + '\t' + end + '\t' + feature + '\t' + reads + '\t' +
                                     tas + '\n')
