conservation = read.table('~/EnTrI/results/homologs.count', row.names = 1)
sorted_conservation <- conservation[ order(row.names(conservation)), ]
clusts <- sort(row.names(conservation))
essentiality = read.table('~/EnTrI/results/homologs.essentiality', row.names = 1)
sorted_essentiality <- essentiality[ order(row.names(essentiality)), ]
plot(sorted_conservation,sorted_essentiality, pch=20, xlab='Conservation level', ylab = 'Average insertion index')
# lines(loess.smooth(sorted_conservation,sorted_essentiality, span=0.2), col=2, lwd=5)
# abline(lm(sorted_essentiality~sorted_conservation), col='red', lwd=3)
