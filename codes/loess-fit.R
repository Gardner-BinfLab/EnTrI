indir <- "../results/per-species-ii/"
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names
list_of_files <- list.files(path=indir, full.names=T, recursive=FALSE)
for (filename in list_of_files)
{
  iitable = read.table(filename, header = FALSE)
  ii = iitable[,2][iitable[,2]>=0]
  h <- hist(ii, breaks =seq(min(ii),max(ii)+1,0.02), plot=FALSE)
  plot(h)
  l = loess(h$counts ~ h$mids)
  lines(l, col="red")
}