from re import match

old_annot = '../data/embl/chromosome_old/Enterobacter_cloacae_subsp_cloacae_NCTC_9394_v1.embl'
prokka_annot = '../results/prokka/ENC/ENC.embl'
new_annot = '../data/embl/chromosome/Enterobacter_cloacae_subsp_cloacae_NCTC_9394_v1.embl'

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
                    locus_tag = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
                    flag = 2
            else:
                break

oldline = ''
addline = ''
counter = 1
with open(old_annot, 'r') as oldfile:
    with open(prokka_annot, 'r') as addfile:
        with open(new_annot, 'a') as newfile:
            while not oldline.startswith('SQ') and not addline.startswith('SQ'):
                while not match('FT\s+CDS\s+', oldline) and not oldline.startswith('SQ'):
                    oldline = oldfile.readline()
                while not match('FT\s+CDS\s+', addline) and not addline.startswith('SQ'):
                    addline = addfile.readline()
                if oldline.startswith('SQ') or addline.startswith('SQ'):
                    break
                start_end = match('FT\s+CDS\s+[^\d]*(\d+)\.\.>?(\d+)', oldline)
                startold = int(start_end.group(1))
                endold = int(start_end.group(2))
                start_end = match('FT\s+CDS\s+[^\d]*(\d+)\.\.(\d+)', addline)
                startadd = int(start_end.group(1))
                endadd = int(start_end.group(2))
                # overlap = (min(endadd, endold) - max(endadd, endold)) /
                if startold == startadd or endold == endadd:
                    oldline = oldfile.readline()
                    addline = addfile.readline()
                elif startold < startadd:
                    oldline = oldfile.readline()
                else:
                    newfile.write(addline)
                    addline = addfile.readline()
                    while not match('FT\s+(CDS|tRNA|misc_feature|sig_peptide|assembly_gap|rRNA)\s+', addline):
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