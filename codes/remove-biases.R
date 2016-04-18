library(stringr)
library(gdata)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names

biasespath <- c("../results/check-biases/with-ends.txt", "../results/check-biases/without-ends.txt")
for (item in biasespath)
{
  biasestable = read.table(item, header = FALSE)
  colnames(biasestable) <- c("name", "ii", "dist", "gc")
  j = 1
  prevj = 1
  for (i in seq(1, length(names)-1))
  {
    ii <- c()
    d <- c()
    while((i == length(names) & j <= length(biasestable$name)) | !startsWith(toString(biasestable$name[j]), names[i+1]))
    {
      ii <- c(ii, biasestable$ii[j])
      d <- c(d, biasestable$dist[j])
      j = j + 1
    }
    fit <- gam(ii~s(d))
    gamprediction <- predict.gam(fit, data.frame(d))
    ii = (ii * mean(ii)) / gamprediction
    biasestable$ii[prevj : (j-1)] = ii
    prevj = j
  }
  ii <- c(biasestable$ii)
  gc = c(biasestable$gc)
  fit <- gam(ii~s(gc))
  gamprediction <- predict.gam(fit, data.frame(gc))
  biasestable$ii = (biasestable$ii * mean(biasestable$ii))/gamprediction
  write.table(biasestable, file = paste("../results/remove-biases/normalised-iis-", strsplit(basename(item), "\\.")[[1]][1], ".txt",sep=''), quote = FALSE, row.names = FALSE, col.names = FALSE) 
}

# hist(biasestable$ii, xlim=c(0,3),breaks=1e5)
# http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best