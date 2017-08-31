seq = 'ctgtctcttatacacatctcttctatagtgtcacctaaatagggataacagggtaatgaattcgctcaaaatctctgatgttacattgcacaagataaaaatatatcatcatgaacaataaaactgtctgcttacataaacagtaatacaaggggtgttatgagccatattcaacgggaaacgtcttgctcgaggccgcgattaaattccaacatggatgctgatttatatgggtataaatgggctcgcgataatgtcgggcaatcaggtgcgacaatctatcgattgtatgggaagcccgatgcgccagagttgtttctgaaacatggcaaaggtagcgttgccaatgatgttacagatgagatggtcagactaaactggctgacggaatttatgcctcttccgaccatcaagcattttatccgtactcctgatgatgcatggttactcaccactgcgatccccggaaaaacagcattccaggtattagaagaatatcctgattcaggtgaaaatattgttgatgcgctggcagtgttcctgcgccggttgcattcgattcctgtttgtaattgtccttttaacagcgatcgcgtatttcgtctcgctcaggcgcaatcacgaatgaataacggtttggttgatgcgagtgattttgatgacgagcgtaatggctggcctgttgaacaagtctggaaagaaatgcataaacttttgccattctcaccggattcagtcgtcactcatggtgatttctcacttgataaccttatttttgacgaggggaaattaataggttgtattgatgttggacgagtcggaatcgcagaccgataccaggatcttgccatcctatggaactgcctcggtgagttttctccttcattacagaaacggctttttcaaaaatatggtattgataatcctgatatgaataaattgcagtttcatttgatgctcgatgagtttttctaatcagaattggttaattggttgtaacactggcagagcattacgctgacttgacgggacggcggctttgttgaataaatcgaacttttgctgagttgaaggatcagatcacgcatcttcccgacaacgcagaccgttccgtggcaaagcaaaagttcaaaatcaccaactggtccacctacaacaaagctctcatcaaccgtggcggggatcctctagagtcgactggcaaacagctattatgggtattatgggtaatacgactcactatagggagatgtgtataagagacag'
codon = {'a':'t','t':'a','g':'c','c':'g'}
revcom = ''
for item in seq:
    revcom +=codon[item]

revcom=revcom[::-1]

index_stop_endfor = [-1,-1,-1]
index_start_endfor = [-1,-1,-1]
for i in range(0,len(seq),3):
    for j in range(0,3):
        if seq[i+j:i+j+3] == 'tag' or seq[i+j:i+j+3] == 'taa' or seq[i+j:i+j+3] == 'tga':
            index_stop_endfor[j] = i+j
        if seq[i+j:i+j+3] == 'atg':
            index_start_endfor[j] = i+j

print('\nForward end:')
for i in range(0,3):
    if index_start_endfor[i] > index_stop_endfor[i]:
        print(str(i) + '\t' + str(index_start_endfor[i]))
        print(seq)

# index_stop_begfor = [-1,-1,-1]
# index_start_begfor = [-1,-1,-1]
# for i in range(len(seq)-3,-1,-3):
#     for j in range(0,3):
#         if seq[i+j:i+j+3] == 'tag' or seq[i+j:i+j+3] == 'taa' or seq[i+j:i+j+3] == 'tga':
#             index_stop_begfor[j] = i+j
#         if seq[i+j:i+j+3] == 'atg':
#             index_start_begfor[j] = i+j
#
# print('\nForward start:')
# for i in range(0,3):
#     if index_start_begfor[i] < index_stop_begfor[i]:
#         print('+' + str(i) + '\t')

# index_stop_begrev = [-1,-1,-1]
# index_start_begrev = [-1,-1,-1]
# for i in range(len(revcom)-3,-1,-3):
#     for j in range(0,3):
#         if revcom[i+j:i+j+3] == 'tag' or revcom[i+j:i+j+3] == 'taa' or revcom[i+j:i+j+3] == 'tga':
#             index_stop_begrev[j] = i+j
#         if revcom[i+j:i+j+3] == 'atg':
#             index_start_begrev[j] = i+j
#
# print('\nReverse start:')
# for i in range(0,3):
#     if index_start_begrev[i] < index_stop_begrev[i]:
#         print('-' + str(i) + '\t')

index_stop_endrev = [-1,-1,-1]
index_start_endrev = [-1,-1,-1]
for i in range(0,len(revcom),3):
    for j in range(0,3):
        if revcom[i+j:i+j+3] == 'tag' or revcom[i+j:i+j+3] == 'taa' or revcom[i+j:i+j+3] == 'tga':
            index_stop_endrev[j] = i+j
        if revcom[i+j:i+j+3] == 'atg':
            index_start_endrev[j] = i+j

print('\nReverse end:')
for i in range(0,3):
    if index_start_endrev[i] > index_stop_endrev[i]:
        print(str(i) + '\t' + str(index_start_endrev[i]))
        print(revcom)

print('\nlength: ' + str(len(seq)))