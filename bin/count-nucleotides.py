from os import listdir

genomedir = '../data/fasta-genome/chromosome/'
list_of_files = listdir(genomedir)
list_of_files.remove('b.fa')
g = 0
c = 0
a = 0
t = 0
for filename in list_of_files:
    with open(genomedir+filename) as fromfile:
        line = fromfile.readline()
        line = fromfile.readline()
        while line and not line.startswith('>'):
            g += line.count('g')
            c += line.count('c')
            a += line.count('a')
            t += line.count('t')
            line = fromfile.readline()
total = a + g + c + t
ap = a * 100 / total
gp = g * 100 / total
cp = c * 100 / total
tp = t * 100 / total
print('G% =  ' + str(gp))
print('C% =  ' + str(cp))
print('A% =  ' + str(ap))
print('T% =  ' + str(tp))