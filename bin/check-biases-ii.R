library(stringr)
library("MASS")
library(gdata)
library(mgcv)
library("dbscan")
library("ROCR")
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261", "BW25113", "EC958")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi", "Salmonella typhimurium SL3261",
         "Escherichia coli BW25113", "Escherichia coli ST131 EC958")
names(dict) <- names
sp = 0.2
biasespath <- "../results/biases/check-biases/"
outdir <- "../results/biases/dbscan"
list_of_files <- list.files(path=biasespath, full.names=T, recursive=FALSE)
iitotal = c()
gctotal = c()
dtotal = c()
ii_dnormalisedtotal = c()
ii_dgcnormalisedtotal = c()
pdf("../figures/biases.pdf")
for (filename in list_of_files)
{
  biasestable = read.table(filename, header = FALSE, stringsAsFactors = FALSE)
  colnames(biasestable) <- c("name", "ii", "essentiality", "dist", "gc", "length")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  ii <- c(biasestable$ii)
  iitotal <- c(iitotal, ii)
  gc <- c(biasestable$gc)
  gctotal <- c(gctotal , gc)
  d <- c(biasestable$dist)
  dtotal <- c(dtotal , d)
  
  plot(d, ii, pch = '.', ylim=c(0,5), xlab = "Gene position", ylab = "insertion index",
       main = paste("Distance bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2)
  # abline(mean(ii),0,col='green', lwd=5)
  lines(loess.smooth(d,ii, span=sp), col=2, lwd=5)
  
  plot(gc, ii, pch = '.', ylim=c(0,5), xlab = "GC content", ylab = "insertion index",
       main = paste("GC bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2)
  # abline(mean(ii),0,col='green', lwd=5)
  lines(loess.smooth(gc,ii, span=sp), col=2, lwd=5)
  
  fit <- loess(ii~d, span = sp)
  loessprediction <- predict(fit, d)
  ii_dnormalised = ii/(loessprediction/mean(ii))
  ii_dnormalisedtotal <- c(ii_dnormalisedtotal, ii_dnormalised)
  plot(d, ii_dnormalised, pch='.', xlab = "Gene position", ylab = "insertion index",
       main = paste(dict[strsplit(basename(filename), "\\.")[[1]][1]], '- normalised distance'), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,5))
  lines(loess.smooth(d,ii_dnormalised, span=sp), col=2, lwd=5)
  
  fit <- loess(ii_dnormalised~gc, span =sp)
  loessprediction <- predict(fit, gc)
  ii_dgcnormalised = ii_dnormalised/(loessprediction/mean(ii_dnormalised))
  ii_dgcnormalisedtotal <- c(ii_dgcnormalisedtotal, ii_dnormalised)
  plot(gc, ii_dgcnormalised, pch='.', xlab = "GC content", ylab = "insertion index",
       main = paste(dict[strsplit(basename(filename), "\\.")[[1]][1]], '- normalised GC'), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,5))
  lines(loess.smooth(gc,ii_dgcnormalised, span=sp), col=2, lwd=5)
}

plot(gctotal, iitotal, pch = '.', ylim=c(0,5), xlab = "GC content", ylab = "insertion index", main = "GC bias", 
     cex.lab = 2, cex.axis = 2, cex.main =2)
lines(loess.smooth(gctotal,iitotal, span=sp), col=2, lwd=5)

dev.off()

pdf("~/EnTrI/results/insertion-indices/results-normalised.pdf")
s = 1
for (filename in list_of_files)
{
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  biasestable = read.table(filename, header = FALSE, stringsAsFactors = FALSE)
  colnames(biasestable) <- c("name", "ii", "essentiality", "dist", "gc")
  biasestable$ii <- ii_dgcnormalisedtotal[s:(s+length(biasestable$ii)-1)]
  # for (i in seq(1,length(biasestable$ii)))
  # {
  #   if (biasestable$ii[i] < 0)
  #   {
  #     print(biasestable$ii[i])
  #     biasestable$ii[i] = 0
  #   }
  # }
  ii=biasestable$ii
  # nG = length(ii)
  # 
  # #identify second maxima
  # h <- hist(ii, breaks=0:(max(ii)*50+1)/50,plot=FALSE)
  # maxindex <- which.max(h$density[3:length(h$density)])
  # maxval <- h$mids[maxindex+2]
  # 
  # #find inter-mode minima with loess
  # r <- floor(maxval *1000)
  # I = ii < r / 1000
  # h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
  # lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  # m1 = h1$mids[which.min(predict(lo))]
  # m2 = h$mids[max(which(h$counts>5))]
  # I1 = ((ii < m1)&(ii > 0))
  # I2 = ((ii >= m1)&(ii < m2))
  # 
  # f1 = (sum(I1) + sum(ii == 0))/nG
  # f2 = (sum(I2))/nG
  # 
  # # if (strsplit(basename(filename),"\\.")[[1]][1] == 'SEN')
  # # {
  # #   m3 = h$mids[min(which(diff(sign(diff(h$counts)))==-2)+1)] # the index of the first local maximum
  # #   I1 = ((ii < m1)&(ii >= m3))
  # # }
  # # 
  # d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
  # d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves
  # 
  # #plots
  # if (strsplit(basename(filename),"\\.")[[1]][1] != 'SEN')
  # {
  #   hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  #   lines(0:2000/500, f1*dgamma(0:2000/500, 1, d1$estimate[2])) # was [2]
  #   lines(0:2000/500, f2*dgamma(0:2000/500, d2$estimate[1], d2$estimate[2]))
  #   lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
  #   upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
  #   essen <- lower/10000
  #   ambig <- upper/10000
  #   noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  # }
  # else
  # {
  #   hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  #   lines(0:2000/500, f1*dgamma(0:2000/500, d1$estimate[1], d1$estimate[2])) # was [2]
  #   lines(0:2000/500, f2*dgamma(0:2000/500, d2$estimate[1], d2$estimate[2]))
  #   lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
  #   upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
  #   essen <- lower/10000
  #   ambig <- upper/10000
  #   noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  # }
  # 
  # lines(c(essen, essen), c(0,20), col="red")
  # lines(c(ambig, ambig), c(0,20), col="red")
  # lines(c(noness, noness), c(0,20), col="red")
  # 
  # mtext(paste(essen, ":", "Essential changepoint"), side=3, adj=1, padj=2)
  # mtext(paste(ambig, ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
  # mtext(paste(noness, ":", "Non-essential changepoint"), side=3, adj=1, padj=5.5)
  # 
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  res <- dbscan(as.matrix(ii), minPts = 200, eps = 0.05)
  # res <- optics(as.matrix(ii), minPts = 200, eps = 0.05)
  # res <- extractDBSCAN(res, 0.05)
  plot(ii,col=res$cluster+1)
  ess <- res$cluster[which.min(ii)]
  print(sum(res$cluster == ess))
  nes <- getmode(res$cluster)
  belthr <- max(ii[res$cluster==nes])
  essthr <- max(ii[res$cluster==ess])
  nesthr <- min(ii[res$cluster==nes])
  
  
  for (i in (1:length(ii)))
  {
    if (res$cluster[i] == ess)
    {
      biasestable$essentiality[i] = "essential"
    } else if(res$cluster[i] == nes)
    {
      biasestable$essentiality[i] = "non-essential"
    } else if(ii[i] > belthr)
    {
      biasestable$essentiality[i] = "beneficial-loss"
    } else
    {
      biasestable$essentiality[i] = "ambiguous"
    }
  }
  
  hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  lines(c(essthr, essthr), c(0,20), col="red")
  lines(c(nesthr, nesthr), c(0,20), col="red")
  lines(c(belthr, belthr), c(0,20), col="red")
  mtext(paste(round(essthr,digits = 3), ":", "Essential changepoint"), side=3, adj=1, padj=2)
  mtext(paste(round(nesthr,digits = 3), ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
  mtext(paste(round(belthr,digits = 3), ":", "Non-essential changepoint"), side=3, adj=1, padj=5.5)
  
  esslevel <- ifelse(ii > 0, log(ii/essthr), -4.5)
  
  write.table(cbind(biasestable[1:3], esslevel), file = paste(outdir, '/', basename(filename), sep=''), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
  s = s + length(biasestable$ii)
}

dev.off()