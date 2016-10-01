library(stringr)
library("MASS")
require(reshape2)  # this is the library that lets you flatten out data
require(ggplot2)
iis_path <- "../results/define-core-accessory-hieranoid-ancestralii/R/iis.dat"
essenpath <- "../results/define-core-accessory-hieranoid-ancestralii/R/essentiality-threshold.dat"


ii = as.double(unlist(read.table(iis_path)))
nG = length(ii)

#identify second maxima
h <- hist(ii, breaks=0:(max(ii)*50+1)/50,plot=FALSE)
maxindex <- which.max(h$density[3:length(h$density)])
maxval <- h$mids[maxindex+2]

#find inter-mode minima with loess
r <- floor(maxval *1000)
I = ii < r / 1000
h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density

m1 = h1$mids[which.min(predict(lo))]
m2 = h$mids[max(which(h$counts>5))]
I1 = ((ii < m1)&(ii > 0))
I2 = ((ii >= m1)&(ii < m2))

f1 = (sum(I1) + sum(ii == 0))/nG
f2 = (sum(I2))/nG

d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves
h <- hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", plot=FALSE)
lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
essen <- lower/10000

write.table(essen, file=essenpath, quote = FALSE, col.names = FALSE, row.names = FALSE)