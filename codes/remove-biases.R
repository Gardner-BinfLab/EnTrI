library(stringr)
library(gdata)
biasespath <- "../results/check-biases.out"
biasestable = read.table(biasespath, header = FALSE)
colnames(biasestable) <- c("name", "ii", "dist", "gc")

names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names
j = 1
prevj = 1
for (i in seq(1, length(names)-1))
{
  ii <- c()
  d <- c()
  while((i == length(names) & j <= length(biasestable$name)) | !startsWith(biasestable$name[j], names[i+1]))
  {
    ii <- c(ii, biasestable$ii[j])
    d <- c(d, biasestable$dist[j])
    j = j + 1
  }
  l = loess(ii ~ d)
  loessprediction <- predict(l, d)
  minimum = min(loessprediction[loessprediction>0])
  for (k in seq(1, length(loessprediction)))
  {
    if (loessprediction[k] <= 0)
    {
      biasestable$ii[k] = minimum + biasestable$ii[k] - loessprediction[k]
      loessprediction[k] = minimum
    }
  }
  ii = (ii * mean(ii)) / loessprediction
  biasestable$ii[prevj : (j-1)] = ii
  prevj = j
}
l = loess(biasestable$ii ~ biasestable$gc)
loessprediction <- predict(l, biasestable$gc)
minimum = min(loessprediction[loessprediction>0])
for (k in seq(1, length(loessprediction)))
{
  if (loessprediction[k] <= 0)
  {
    biasestable$ii[k] = minimum + biasestable$ii[k] - loessprediction[k]
    loessprediction[k] = minimum
  }
}
biasestable$ii = (biasestable$ii * mean(biasestable$ii))/loessprediction
write.table(biasestable, file = "../results/normalised-iis.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# hist(biasestable$ii, xlim=c(0,3),breaks=1e5)
# http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best