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
  
  pdf(paste("../results/biases-", strsplit(basename(item), "\\.")[[1]][1], ".pdf",sep=''))
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  
  plot(biasestable$dist, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "Gene position", ylab = "insertion index",
       main = "Distance bias", cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(biasestable$dist, biasestable$ii, span = 0.2), col = 2, lwd=5)
  
  j = 1
  for (i in seq(1, length(names)-1))
  {
    ii <- c()
    d <- c()
    gc <- c()
    while(!startsWith(biasestable$name[j], names[i+1]))
    {
      ii <- c(ii, biasestable$ii[j])
      d <- c(d, biasestable$dist[j])
      gc <- c(gc, biasestable$gc[j])
      j = j + 1
    }
    plot(d, ii, pch = '.', ylim=c(0,3), xlab = "Gene position", ylab = "insertion index",
         main = dict[names[i]], cex.lab = 2, cex.axis = 2, cex.main =2)
    lines(loess.smooth(d, ii, span = 0.2), col = 2, lwd=5)
    
    # normalization by loess and mean
    # abline(mean(ii),0, col=3, lwd=5)
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
    plot(d, ii/(loessprediction/mean(ii)), pch='.', xlab = "Gene position", ylab = "insertion index",
         main = paste(dict[names[i]], '-normalised', sep=''), cex.lab = 2, cex.axis = 2, cex.main =2,
         ylim=c(0,2))
    lines(loess.smooth(d, ii/(loessprediction/mean(ii)), span = 0.2), col = 2, lwd=5)
    
    #   #relation between gc content and distance from origin
    #   plot(d,gc,pch='.')
    #   lines(loess.smooth(d, gc, span = 0.2), col = 2, lwd=5)
    
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
  lines(loess.smooth(d, ii, span = 0.2), col = 2, lwd=5)
  # abline(mean(ii),0, col=3, lwd=5)
  l = loess(ii ~ d, span=0.2)
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
  plot(d, ii/(loessprediction/mean(ii)), pch='.', xlab = "Gene position", ylab = "insertion index",
       main = paste(dict[names[i]], '-normalised', sep=''), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,2))
  lines(loess.smooth(d, ii/(loessprediction/mean(ii)), span = 0.2), col = 2, lwd=5)
  
  plot(biasestable$gc, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "GC content", ylab = "insertion index", main = "GC bias",
       cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(biasestable$gc, biasestable$ii, span = 0.2), col = 2, lwd=5)
  
  # abline(mean(biasestable$ii),0, col=3, lwd=5)
  
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
  plot(biasestable$gc, biasestable$ii/(loessprediction/mean(biasestable$ii)), pch='.', ylim=c(0,3), xlab = "GC content",
       ylab = "insertion index", main = "GC bias-normalised", cex.lab = 2, cex.axis = 2, cex.main =2)
  lines(loess.smooth(biasestable$gc, biasestable$ii/(loessprediction/mean(biasestable$ii)), span = 0.2), col = 2, lwd=5)
  
  dev.off()
}