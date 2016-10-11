library(stringr)
library("MASS")
# args <- commandArgs(trailingOnly = TRUE)
# clusters <- args[1]
clusters_path <- "../results/merge-clust-plot"

list_of_files <- list.files(path=clusters_path, full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261")
genuses = c(1, 2, 3, 2, 2, 4, 4, 5, 5, 5, 5, 5, 5)
names(genuses) <- names
numspecies = length(names)
file_II = list()
file_size = list()
file_group = list()
for (filename in list_of_files)
{
  clustspecies = c()
  cluster <- as.matrix(read.table(filename))
  i_sum = 0
  l_sum = 0
  cluster_size = nrow(cluster)
  clust_with_ii_size = 0
  for (i in (1:cluster_size))
  {
    if (as.numeric(cluster[i,5]) >= 0)
    {
      match = str_match(cluster[i,2], "([[:graph:]]+)\\_[[:alnum:]]+")[2]
      if (is.na(match))
      {
        match = str_match(cluster[i,2], "([[:alpha:]]+)[[:digit:]]+")[2]
      }
      if (match %in% names)
      {
        clustspecies = c(clustspecies, match) 
      }
      i_sum = i_sum + as.numeric(cluster[i, 5])
      l_sum = l_sum + (as.numeric(cluster[i, 4]) - as.numeric(cluster[i, 3]) + 1) * 3
      clust_with_ii_size = clust_with_ii_size  + 1
    }
  }
  if (l_sum > clust_with_ii_size * 60)
  {
    file_II[basename(filename)] = i_sum / clust_with_ii_size
    file_size[basename(filename)] = cluster_size
    unique_clustspecies = unique(clustspecies)
    one_or_more = length(unique_clustspecies)
    greater_than_one = length(table(clustspecies)[table(clustspecies)>1])
    cluster_genuses = c()
    for (item in unique_clustspecies)
    {
      cluster_genuses <- c(cluster_genuses, genuses[unique_clustspecies])
    }
    # if (as.numeric(one_or_more) <= 0.3 * numspecies)
    if (length(unique(cluster_genuses)) <= 1)
    {
      file_group[basename(filename)] = 'ORFan'
    }
    else if (as.numeric(one_or_more) - as.numeric(greater_than_one) >= 0.7 * as.numeric(one_or_more))
    {
      file_group[basename(filename)] = 'Single-copy'
    }
    else
    {
      file_group[basename(filename)] = 'Multiple-copy'
    }
  }
}
insertion_index <- sapply(file_II, function(x){as.numeric(x[1])})
size_index <- sapply(file_size, function(x){as.numeric(x[1])})
group_index <- file_group

pdf("../figures/cluster-essentiality2.pdf")

m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
temp <- split.screen(m)

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

screen(1)
ii = insertion_index
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
ambig <- upper/10000
noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])

cuts <- cut(h$breaks, c(-Inf,essen, ambig, noness, Inf))
plot(h, col=c("darkgoldenrod4", "gray", "turquoise4", "darkmagenta")[cuts], xlab = "Insertion Index", main ="All clusters", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,200,50), labels=c(0,NA,NA,NA,200), cex.axis=1.5)
text(2.85,190, paste("n =", length(ii)), lty=1, lwd=4, cex=1.15, bty="n")
legend(2.4,180, c("Essential","Ambiguous","Non-essential", "Beneficial loss"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.15,
       col=c("darkgoldenrod4","gray","turquoise4", "darkmagenta"), bty="n")


orfans = c()
orfans_non = c()
orfans_es = c()
orfans_ben = c()
single_occurrence = c()
single_non = c()
single_es = c()
single_ben = c()
multiple_copies = c()
multiple_non = c()
multiple_es = c()
multiple_ben = c()
for (item in names(size_index))
{
  if (group_index[item] == 'ORFan')
  {
    orfans = c(orfans, insertion_index[item])
    if (insertion_index[item] < 0.2)
      orfans_es = c(orfans_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      orfans_non = c(orfans_non, insertion_index[item])
    else
      orfans_ben = c(orfans_ben, insertion_index[item])
  } 
  else if (group_index[item] == 'Single-copy')
  {
    single_occurrence = c(single_occurrence, insertion_index[item])
    if (insertion_index[item] < 0.2)
      single_es = c(single_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      single_non = c(single_non, insertion_index[item])
    else
      single_ben = c(single_ben, insertion_index[item])
  }
  else
  {
    multiple_copies = c(multiple_copies, insertion_index[item])
    if (insertion_index[item] < 0.2)
      multiple_es = c(multiple_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      multiple_non = c(multiple_non, insertion_index[item])
    else
      multiple_ben = c(multiple_ben, insertion_index[item])
  }
}

screen(2)
par(mar=c(5.1,2.5,4.1,1))
ii = orfans
h <- hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,essen, ambig, noness, Inf))
plot(h, col=c("darkgoldenrod4", "gray", "turquoise4", "darkmagenta")[cuts], xlab = "", main ="Genus specific", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,200,50), labels=c(0,NA,NA,NA,200), cex.axis=1.5)
text(2.85,190, paste("n =", length(ii)), lty=1, lwd=4, cex=1.15, bty="n")

screen(3)
par(mar=c(5.1,1,4.1,1))
ii = single_occurrence
h <- hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,essen, ambig, noness, Inf))
plot(h, col=c("darkgoldenrod4", "gray", "turquoise4", "darkmagenta")[cuts], xlab = "", main ="Single copy", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,200,50), labels=c(NA,NA,NA,NA,NA), cex.axis=1.5)
text(2.85,190, paste("n =", length(ii)), lty=1, lwd=4, cex=1.15, bty="n")

screen(4)
#par(mar=c(2,1,2,1))
par(mar=c(5.1,1,4.1,1))
ii = multiple_copies
h <- hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,essen, ambig, noness, Inf))
plot(h, col=c("darkgoldenrod4", "gray", "turquoise4", "darkmagenta")[cuts], xlab = "", main ="Multi-copy", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,200,50), labels=c(NA,NA,NA,NA,NA), cex.axis=1.5)
text(2.85,190, paste("n =", length(ii)), lty=1, lwd=4, cex=1.15, bty="n")



close.screen(all.screens = TRUE)

dev.off()