library(stringr)
library(gdata)
biasespath <- "../results/check-biases.out"
biasestable = read.table(biasespath, header = FALSE)
colnames(biasestable) <- c("name", "ii", "dist", "gc")

names = c("ROD", "CS17", "ENC", "Ecoli9000q", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli 9000", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names

pdf("../results/biases.pdf")

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

plot(biasestable$dist, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "Distance from the origin", ylab = "insertion index",
     main = "Distance bias", cex.lab = 2, cex.axis = 2, cex.main =2)
lines(lowess(biasestable$dist, biasestable$ii, f = 0.2), col = 2, lwd=5)

j = 1
for (i in seq(1, length(names)-1))
{
  ii <- c()
  d <- c()
  while(!startsWith(biasestable$name[j], names[i+1]))
  {
    ii <- c(ii, biasestable$ii[j])
    d <- c(d, biasestable$dist[j])
    j = j + 1
  }
  plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Distance from the origin", ylab = "insertion index",
       main = dict[names[i]], cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(lowess(d, ii, f = 0.2), col = 2, lwd=5)
}
ii <- c()
d <- c()
i = i + 1
while(j <= length(biasestable$name))
{
  ii <- c(ii, biasestable$ii[j])
  d <- c(d, biasestable$dist[j])
  j = j + 1
}
plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Distance from the origin", ylab = "insertion index",
     main = dict[names[i]], cex.lab = 2, cex.axis = 2, cex.main =2)
lines(lowess(d, ii, f = 0.2), col = 2, lwd=5)

plot(biasestable$gc, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "GC content", ylab = "insertion index", main = "GC bias",
     cex.lab = 2, cex.axis = 2, cex.main =2)
regln <- lm(biasestable$gc ~ biasestable$ii)
abline(regln, col=3, lwd=5)

dev.off()