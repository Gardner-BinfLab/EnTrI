#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from re import match, findall
from argparse import ArgumentParser
import string
import random


def read_sequence(seqfile):
    seq = ''
    with open(seqfile, 'r') as emblfile:
        emblline = emblfile.readline()
        while match('SQ\s+Sequence', emblline) is None:
            emblline = emblfile.readline()
        while emblline:
            emblline = emblfile.readline()
            spaceless = emblline.replace(' ', '')
            match_res = match('([A-Za-z]*)\d+', spaceless)
            if match_res:
                digitless = match_res.group(1)
                lower_case = digitless.lower()
                seq += lower_case.replace('u', 't')
    return seq


def translate(gene, comp):
    translation = ''
    if comp == 'Complement':
        gene = gene[::-1]
        gene = gene.replace('a', 'p')
        gene = gene.replace('t', 'q')
        gene = gene.replace('c', 'r')
        gene = gene.replace('g', 's')
        gene = gene.replace('p', 't')
        gene = gene.replace('q', 'a')
        gene = gene.replace('r', 'g')
        gene = gene.replace('s', 'c')
    for j in range(0, int(len(gene))-3):
        if not(j % 60):
            translation += '\n'
        translation += gene[j]
    j += 1
    if gene[j:j+3] != 'TAG' and gene[j:j+3] != 'TAA' and gene[j:j+3] != 'TGA':
        if not(j % 60):
            translation += '\n'
        translation += gene[j:j+3]
    return translation

parser = ArgumentParser(description='Converts embl files to DNA fasta.')
parser.add_argument('embl', help='Embl file name (source)')
parser.add_argument('fasta', help='Fasta file name (destination)')
args = parser.parse_args()
source = args.embl
destination = args.fasta

#source = '../protein-clustering/sequences/embl/AE016877.txt'
#destination = '../protein-clustering/sequences/dna-fasta/AE016877.txt'

sequence = read_sequence(source)
t_dict = {}
s, e = [], []
p, nt, t, g, comment, se, tr, ec = '', '', '', '', '', '', '', ''
locus_tag_counter = 0
ID = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase) for _ in range(5))
with open(source, 'r') as from_file:
    with open(destination, 'w') as to_file:
        line = from_file.readline()
        while line:
            if match('ID\s+(\S+);', line):
                match_result = match('ID\s+(\S+);', line)
                ID = match_result.group(1)
            elif match('FT\s+CDS\s+>?<?(\d+)\.\.>?<?(\d+)', line):
                match_result = match('FT\s+CDS\s+>?<?(\d+)\.\.>?<?(\d+)', line)
                if se and 'misc_RNA' not in t:
                    if not tr:
                        seq_to_translate = ''
                        for i in range(0, len(s)):
                            seq_to_translate += sequence[int(s[i])-1:int(e[i])]
                        tr = translate(seq_to_translate, comment)
                    if not t:
                        t = ID[:3] + 'LocTag_' + str(locus_tag_counter)
                        locus_tag_counter += 1
                    to_file.write('>{0} [{1}/{2} ({3})] {4} {5} {6} [EC_number={7}]'.format(t, ID, se, comment, p, nt, g, ec))
                    to_file.write(tr)
                    to_file.write('\n')
                s, e = [], []
                se = ''
                s.append(str(min(int(match_result.group(2)), int(match_result.group(1)))))
                e.append(str(max(int(match_result.group(2)), int(match_result.group(1)))))
                se += s[0] + '-' + e[0]
                comment = 'Forward'
                p, nt, t, g, tr, ec = '', '', '', '', '', ''
            elif match('FT\s+CDS\s+complement\(>?<?(\d+)\.\.>?<?(\d+)\)', line):
                match_result = match('FT\s+CDS\s+complement\(>?<?(\d+)\.\.>?<?(\d+)\)', line)
                if se and 'misc_RNA' not in t:
                    if not tr:
                        seq_to_translate = ''
                        for i in range(0, len(s)):
                            seq_to_translate += sequence[int(s[i])-1:int(e[i])]
                        tr = translate(seq_to_translate, comment)
                    if not t:
                        t = ID[:3] + 'LocTag_' + str(locus_tag_counter)
                        locus_tag_counter += 1
                    to_file.write('>{0} [{1}/{2} ({3})] {4} {5} {6} [EC_number={7}]'.format(t, ID, se, comment, p, nt, g, ec))
                    to_file.write(tr)
                    to_file.write('\n')
                s, e = [], []
                se = ''
                s.append(str(min(int(match_result.group(2)), int(match_result.group(1)))))
                e.append(str(max(int(match_result.group(2)), int(match_result.group(1)))))
                se += s[0] + '-' + e[0]
                comment = 'Complement'
                p, nt, t, g, tr, ec = '', '', '', '', '', ''
            elif match('FT\s+CDS\s+join\([><\d\.,]+\)', line):
                match_result = match('FT\s+CDS\s+join\([><\d\.,]+\)', line)
                if se and 'misc_RNA' not in t:
                    if not tr:
                        seq_to_translate = ''
                        for i in range(0, len(s)):
                            seq_to_translate += sequence[int(s[i])-1:int(e[i])]
                        tr = translate(seq_to_translate, comment)
                    if not t:
                        t = ID[:3] + 'LocTag_' + str(locus_tag_counter)
                        locus_tag_counter += 1
                    to_file.write('>{0} [{1}/{2} ({3})] {4} {5} {6} [EC_number={7}]'.format(t, ID, se, comment, p, nt, g, ec))
                    to_file.write(tr)
                    to_file.write('\n')
                find_result = findall('>?<?\d+\.\.>?<?\d+', line)
                s, e = [], []
                se = ''
                for item in find_result:
                    start_end = match('>?<?(\d+)\.\.>?<?(\d+)', item)
                    s.append(start_end.group(1))
                    e.append(start_end.group(2))
                    if se != '':
                        se += ' ' + s[len(s)-1] + '-' + e[len(s)-1]
                    else:
                        se += s[len(s)-1] + '-' + e[len(s)-1]
                comment = 'Forward'
                p, nt, t, g, tr, ec = '', '', '', '', '', ''
            elif match('FT\s+CDS\s+complement\(join\([><\d\.,]+\)\)', line):
                match_result = match('FT\s+CDS\s+complement\(join\([><\d\.,]+\)\)', line)
                if se and 'misc_RNA' not in t:
                    if not tr:
                        seq_to_translate = ''
                        for i in range(0, len(s)):
                            seq_to_translate += sequence[int(s[i])-1:int(e[i])]
                        tr = translate(seq_to_translate, comment)
                    if not t:
                        t = ID[:3] + 'LocTag_' + str(locus_tag_counter)
                        locus_tag_counter += 1
                    to_file.write('>{0} [{1}/{2} ({3})] {4} {5} {6} [EC_number={7}]'.format(t, ID, se, comment, p, nt, g, ec))
                    to_file.write(tr)
                    to_file.write('\n')
                find_result = findall('>?<?\d+\.\.>?<?\d+', line)
                s, e = [], []
                se = ''
                for item in find_result:
                    start_end = match('>?<?(\d+)\.\.>?<?(\d+)', item)
                    s.append(start_end.group(1))
                    e.append(start_end.group(2))
                    if se != '':
                        se += ' ' + s[len(s)-1] + '-' + e[len(s)-1]
                    else:
                        se += s[len(s)-1] + '-' + e[len(s)-1]
                comment = 'Complement'
                p, nt, t, g, tr, ec = '', '', '', '', '', ''
            elif match('FT\s+\S+\s+\d+\.\.\d+', line) or match('FT\s+\S+\s+complement\(\d+\.\.\d+', line) or\
                    match('FT\s+\S+\s+join\([\d\.,]+', line) or match('FT\s+\S+\s+complement\(join\([\d\.,]+',
                                                                      line):
                line = from_file.readline()
                while not(match('FT\s+\S+\s+\d+\.\.\d+', line) or match('FT\s+\S+\s+complement\(\d+\.\.\d+', line)
                          or match('FT\s+\S+\s+join\([\d\.,]+', line) or
                          match('FT\s+\S+\s+complement\(join\([\d\.,]+', line) or match('SQ\s+Sequence', line)):
                    line = from_file.readline()
                continue
            elif match('FT\s+/gene=\"(\S+)\"', line):
                match_result = match('FT\s+/gene=\"(\S+)\"', line)
                g = '[{}]'.format(match_result.group(1))
            elif match('FT\s+/product=\"(.*)\"', line):
                match_result = match('FT\s+/product=\"(.*)\"', line)
                p = '[{}]'.format(match_result.group(1))
            elif match('FT\s+/product=\"(.*)', line):
                match_result = match('FT\s+/product=\"(.*)', line)
                p = '[{}'.format(match_result.group(1))
                line = from_file.readline()
                while not match('FT\s+(.*)\"', line):
                    match_result = match('FT\s+(.*)', line)
                    p += ' {}'.format(match_result.group(1))
                    line = from_file.readline()
                match_result = match('FT\s+(.*)\"', line)
                p += ' {}]'.format(match_result.group(1))
            elif match('FT\s+/locus_tag=\"(\S+)\"', line):
                match_result = match('FT\s+/locus_tag=\"(\S+)\"', line)
                t = match_result.group(1)
                if t in t_dict.keys():
                    t_dict[t] += 1
                    t = t + '_' + str(t_dict[t])
                else:
                    t_dict[t] = 1
            elif match('FT\s+/EC_number=\"(\S+)\"', line):
                match_result = match('FT\s+/EC_number=\"(\S+)\"', line)
                if not ec == '':
                    ec += ', '
                ec += match_result.group(1)
            elif match('FT\s+/note=\"(.*)\"', line):
                match_result = match('FT\s+/note=\"(.*)\"', line)
                nt = '[{}]'.format(match_result.group(1))
            elif match('FT\s+/note=\"(.*)', line):
                match_result = match('FT\s+/note=\"(.*)', line)
                nt = '[{}'.format(match_result.group(1))
                line = from_file.readline()
                while not match('FT\s+(.*)\"', line):
                    match_result = match('FT\s+(.*)', line)
                    nt += ' {}'.format(match_result.group(1))
                    line = from_file.readline()
                match_result = match('FT\s+(.*)\"', line)
                nt += ' {}]'.format(match_result.group(1))
            line = from_file.readline()
        if se:
            if not tr:
                seq_to_translate = ''
                for i in range(0, len(s)):
                    seq_to_translate += sequence[int(s[i])-1:int(e[i])]
                tr = translate(seq_to_translate, comment)
            if not t:
                t = ID[:3] + 'LocTag_' + str(locus_tag_counter)
            to_file.write('>{0} [{1}/{2} ({3})] {4} {5} {6} [EC_number={7}]'.format(t, ID, se, comment, p, nt, g, ec))
            to_file.write(tr)
            to_file.write('\n')