both = read.table('../results/insertion-indices/insertion-position-bias-first100.out', row.names = 1)
bothsum=colSums(both)
barplot(bothsum[1:24])

posit = read.table('../results/insertion-indices/insertion-position-bias-first100-positive.out', row.names = 1)
positsum=colSums(posit)
barplot(positsum[1:24])

rever = read.table('../results/insertion-indices/insertion-position-bias-first100-reverse.out', row.names = 1)
reversum=colSums(rever)
barplot(reversum[1:24])