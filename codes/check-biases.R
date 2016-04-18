library(stringr)
library(gdata)
library(mgcv)
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
  
  # pdf(paste("../results/biases-", strsplit(basename(item), "\\.")[[1]][1], ".pdf",sep=''))
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  ii <- c(biasestable$ii)
  d <- c(biasestable$dist)
  plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Gene position", ylab = "insertion index",
       main = "Distance bias", cex.lab = 2, cex.axis = 2, cex.main =2)
  fit <- gam(ii~s(d))
  gamprediction <- predict.gam(fit, data.frame(d))
  lines(loess.smooth(d,ii, span=0.2), col=2, lwd=5)
  
  j = 1
  for (i in seq(1, length(names)-1))
  {
    ii <- c()
    d <- c()
    gc <- c()
    while(!startsWith(toString(biasestable$name[j]), names[i+1]))
    {
      ii <- c(ii, biasestable$ii[j])
      d <- c(d, biasestable$dist[j])
      gc <- c(gc, biasestable$gc[j])
      j = j + 1
    }
    plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Gene position", ylab = "insertion index",
         main = dict[names[i]], cex.lab = 2, cex.axis = 2, cex.main =2)
    fit <- gam(ii~s(d))
    gamprediction <- predict.gam(fit, data.frame(d))
    lines(loess.smooth(d,ii, span=0.2), col=2, lwd=5)
    ii = ii/(gamprediction/mean(ii))
    plot(d, ii, pch='.', xlab = "Gene position", ylab = "insertion index",
         main = paste(dict[names[i]], '-normalised', sep=''), cex.lab = 2, cex.axis = 2, cex.main =2,
         ylim=c(0,2))
    fit <- gam(ii~s(d))
    gamprediction <- predict.gam(fit, data.frame(d))
    lines(loess.smooth(d,ii, span=0.2), col=2, lwd=5)
  }
  ii <- c()
  d <- c()
  gc <- c()
  i = i + 1
  while(j <= length(biasestable$name))
  {
    ii <- c(ii, biasestable$ii[j])
    d <- c(d, biasestable$dist[j])
    gc <- c(gc, biasestable$gc[j])
    j = j + 1
  }
  plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Gene position", ylab = "insertion index",
       main = dict[names[i]], cex.lab = 2, cex.axis = 2, cex.main =2)
  fit <- gam(ii~s(d))
  gamprediction <- predict.gam(fit, data.frame(d))
  lines(loess.smooth(d,ii, span=0.2), col=2, lwd=5)
  ii = ii/(gamprediction/mean(ii))
  plot(d, ii, pch='.', xlab = "Gene position", ylab = "insertion index",
       main = paste(dict[names[i]], '-normalised', sep=''), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,2))
  fit <- gam(ii~s(d))
  gamprediction <- predict.gam(fit, data.frame(d))
  lines(loess.smooth(d,ii, span=0.2), col=2, lwd=5)
  
  gc <- c(biasestable$gc)
  ii <- c(biasestable$ii)
  plot(gc, ii, pch = '.', ylim=c(0,3), xlab = "GC content", ylab = "insertion index", main = "GC bias",
       cex.lab = 2, cex.axis = 2, cex.main =2)
  fit <- gam(ii~s(gc))
  gamprediction <- predict.gam(fit, data.frame(gc))
  lines(loess.smooth(gc,ii, span=0.2), col=2, lwd=5)
  ii = ii/(gamprediction/mean(ii))
  plot(gc, ii, pch='.', ylim=c(0,3), xlab = "GC content",
       ylab = "insertion index", main = "GC bias-normalised", cex.lab = 2, cex.axis = 2, cex.main =2)
  fit <- gam(ii~s(gc))
  gamprediction <- predict.gam(fit, data.frame(gc))
  lines(loess.smooth(gc,ii, span=0.2), col=2, lwd=5)
  
  # dev.off()
}