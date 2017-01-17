library(stringr)
library("MASS")
library(gdata)
library(mgcv)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261", "BW25113", "EC958")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi", "Salmonella typhimurium SL3261",
         "Escherichia coli BW25113", "Escherichia coli ST131 EC958")
names(dict) <- names
sp = 0.2
biasespath <- "../results/biases/check-biases/"
outdir <- "../results/biases/normalised-pca"
dir.create(outdir)
list_of_files <- list.files(path=biasespath, full.names=T, recursive=FALSE)
iitotal = c()
gctotal = c()
dtotal = c()
ltotal = c()
ii_dnormalisedtotal = c()
ii_dgcnormalisedtotal = c()
lengths=c()
pdf("../figures/biases.pdf")
for (filename in list_of_files)
{
  biasestable = read.table(filename, header = FALSE)
  colnames(biasestable) <- c("name", "ii", "essentiality", "dist", "gc", "length")
  for (i in seq(1, nrow(biasestable), 1))
  {
    if (biasestable$"gc"[i] < 0.45)
    {
      lengths[as.character(biasestable$"name"[i])] = biasestable$"length"[i]
    }
  }
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  ii <- c(biasestable$ii)
  iitotal <- c(iitotal, ii)
  gc <- c(biasestable$gc)
  gctotal <- c(gctotal , gc)
  d <- c(biasestable$dist)
  dtotal <- c(dtotal , d)
  l <- c(biasestable$length)
  ltotal <- c(ltotal , l)
  
  plot(gc, ii, pch = '.', xlab = "GC content", ylab = "insertion index",
       main = paste("GC bias -", strsplit(basename(filename), "\\.")[[1]][1]), cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(gc,ii, span=sp), col=2, lwd=5)
  
  plot(d, ii, pch = '.', xlab = "Gene position", ylab = "insertion index",
       main = paste("Distance bias -", strsplit(basename(filename), "\\.")[[1]][1]), cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(d,ii, span=sp), col=2, lwd=5)
  
  fit <- loess(ii~d, span=sp)
  loessprediction <- predict(fit, d)
  # ii_dnormalised = ii/(loessprediction/mean(ii))
  ii_dnormalised = ii-(loessprediction-mean(loessprediction))
  
  fit <- loess(ii_dnormalised~gc, span=sp)
  loessprediction <- predict(fit, gc)
  # ii_dgcnormalised = ii_dnormalised/(loessprediction/mean(ii_dnormalised))
  ii_dgcnormalised = ii_dnormalised-(loessprediction-mean(loessprediction))
  
  ii_dnormalisedtotal <- c(ii_dnormalisedtotal, ii_dnormalised)
  ii_dgcnormalisedtotal <- c(ii_dgcnormalisedtotal, ii_dgcnormalised)
  plot(d, ii_dnormalised, pch='.', xlab = "Gene position", ylab = "insertion index",
       main = paste(strsplit(basename(filename), "\\.")[[1]][1], '- normalised distance'), cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(d,ii_dnormalised, span=sp), col=2, lwd=5)
  # fit <- loess(ii_dnormalised~d)
  # loessprediction <- predict(fit, d)
  biasestable$ii = (ii_dgcnormalised- mean(ii_dgcnormalised))/sd(ii_dgcnormalised- mean(ii_dgcnormalised))
  cutoff = 1.644854
  for (i in (1:length(biasestable$ii)))
  {
    if (biasestable$ii[i] > cutoff)
    {
      biasestable$essentiality[i] = "essential"
    }
    else if(biasestable$ii[i] < -cutoff)
    {
      biasestable$essentiality[i] = "beneficial-loss"
    }
    else
    {
      biasestable$essentiality[i] = "non-essential"
    }
  }
  write.table(biasestable[1:3], file = paste(outdir, '/', basename(filename), sep=''), quote = FALSE, col.names = FALSE,
              row.names = FALSE, sep = '\t')
}

plot(gctotal, iitotal, pch = '.', xlab = "GC content", ylab = "insertion index", main = "GC bias", 
     cex.lab = 2, cex.axis = 2, cex.main =2)
lines(loess.smooth(gctotal,iitotal, span=sp), col=2, lwd=5)
# fit <- loess(iitotal~gctotal)
# loessprediction <- predict(fit, gctotal)


fit <- loess(ii_dnormalisedtotal~gctotal)
loessprediction <- predict(fit, gctotal)
iinormalisedtotal = ii_dnormalisedtotal-(loessprediction-mean(loessprediction))
plot(gctotal, iinormalisedtotal, pch='.', xlab = "GC content",
     ylab = "insertion index", main = "GC bias - normalised", cex.lab = 2, cex.axis = 2, cex.main =2)
lines(loess.smooth(gctotal,iinormalisedtotal, span=sp), col=2, lwd=5)
# fit <- loess(iinormalisedtotal~gctotal)
# loessprediction <- predict(fit, gctotal)

plot(ltotal, iitotal, pch = '.', xlab = "Length", xlim=c(0,3000), ylab = "insertion index", main = "Length bias", 
     cex.lab = 2, cex.axis = 2, cex.main =2)
lines(loess.smooth(ltotal,iitotal, span=sp), col=2, lwd=5)

dev.off()