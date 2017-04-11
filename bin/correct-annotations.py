from re import match
from os import listdir

list_of_files = listdir('../data/embl/chromosome_old/')
for genome_name in list_of_files:
    print(genome_name)
    old_annot = '../data/embl/chromosome_old/' + genome_name
    new_annot = '../data/embl/chromosome/' + genome_name
    flag = 0
    with open(old_annot, 'r') as fromfile:
        with open(new_annot, 'w') as tofile:
            for line in fromfile:
                if not line.startswith('SQ'):
                    tofile.write(line)
                    if not flag and match('FT\s+CDS\s+', line):
                        flag = 1
                    if flag == 1 and match('FT\s+/locus_tag=\"(\S+)\"', line):
                        gene = match('FT\s+/locus_tag=\"(\S+)\"', line).group(1)
                        locus_tag = match('([a-zA-Z0-9]+_+|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
                        flag = 2
                else:
                    break

    prokka_annot = '../results/prokka/' + locus_tag + '/' + locus_tag + '.embl'

    oldline = ''
    addline = ''
    nocdsline = ''
    counter = 1
    with open(old_annot, 'r') as oldfile:
        with open(prokka_annot, 'r') as addfile:
            with open(new_annot, 'a') as newfile:
                while not oldline.startswith('SQ') and not addline.startswith('SQ'):
                    while not match('FT\s{3}[a-zA-Z_]+\s+[^\d]*\d+[>|<]?\.\.[>|<]?\d+', oldline) or \
                            match('FT\s{3}source+\s+[^\d]*\d+[>|<]?\.\.[>|<]?\d+', oldline):
                    # while not match('FT\s+CDS\s+[^\d]*\d+[>|<]?\.\.[>|<]?\d+', oldline):
                        oldline = oldfile.readline()
                    while not match('FT\s+CDS\s+[^\d]*(?:FN649414:)?\d+\.\.\d+', addline) and not addline.startswith(
                            'SQ'):
                        addline = addfile.readline()
                    if oldline.startswith('SQ') or addline.startswith('SQ'):
                        break
                    start_end = match(
                        'FT\s{3}[a-zA-Z_]+\s+([^\d]*)(\d+)[>|<]?\.\.[>|<]?(?:\d+,\d+[>|<]?\.\.[>|<]?)*(\d+)', oldline)
                    # start_end = match(
                    #     'FT\s+CDS\s+([^\d]*)(\d+)[>|<]?\.\.[>|<]?(?:\d+,\d+[>|<]?\.\.[>|<]?)*(\d+)', oldline)
                    textold = start_end.group(1)
                    startold = int(start_end.group(2))
                    endold = int(start_end.group(3))
                    if textold.startswith('join('):
                        textold = textold[5:]
                    start_end = match('FT\s{3}CDS\s+([^\d]*)(?:FN649414:)?(\d+)\.\.(\d+)', addline)
                    textadd = start_end.group(1)
                    startadd = int(start_end.group(2))
                    endadd = int(start_end.group(3))
                    if textadd.startswith('join('):
                        textadd = textadd[5:]
                    overlap = (min(endadd, endold) - max(startadd, startold))  # / min(endadd-startadd, endold-startold)
                    if overlap > 0 or endadd - startadd < 300:
                        addline = addfile.readline()
                    elif startold < startadd:
                        oldline = oldfile.readline()
                    else:
                        addline = addline.replace('FN649414:', '')
                        newfile.write(addline)
                        addline = addfile.readline()
                        while not match('FT\s+(CDS|tRNA|misc_feature|sig_peptide|assembly_gap|rRNA)\s+',
                                        addline) and not \
                                addline.startswith('SQ'):
                            match_result = match('(FT\s+/locus_tag=\")\S+\"', addline)
                            if match_result:
                                addline = match_result.group(1)
                                addline = addline + locus_tag + '_' + str(counter) + 'added\"\n'
                                counter += 1
                            newfile.write(addline)
                            addline = addfile.readline()
                while not addline.startswith('SQ'):
                    match_result = match('(FT\s+/locus_tag=\")\S+\"', addline)
                    if match_result:
                        addline = match_result.group(1)
                        addline = addline + locus_tag + '_' + str(counter) + 'added\"\n'
                        counter += 1
                    newfile.write(addline)
                    addline = addfile.readline()

    flag = 0
    with open(old_annot, 'r') as fromfile:
        with open(new_annot, 'a') as tofile:
            for line in fromfile:
                if flag:
                    tofile.write(line)
                elif line.startswith('SQ'):
                    tofile.write(line)
                    flag = 1