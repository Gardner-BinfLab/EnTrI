library(stringr)
list_of_files <- list.files(path="../results/homclust/EFam-clusters", pattern="*.txt", full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "b")
numspecies = length(names)
clusters = matrix(nrow=length(list_of_files), ncol=3)
i = 1
for (filename in list_of_files)
{
  clust_tbl =  read.table(filename, header = FALSE)
  clustspecies = c()
  count = 0
  for (item in clust_tbl[,2])
  {
    count = count + 1
    match = str_match(item, "([[:graph:]]+)\\_[[:alnum:]]+")[2]
    if (is.na(match))
    {
      match = str_match(item, "([[:alpha:]]+)[[:digit:]]+")[2]
    }
    if (match %in% names)
    {
      clustspecies = c(clustspecies, match) 
    }
  }
  clustspecies = unique(clustspecies)
  clusters[i,1]= basename(filename)
  clusters[i,2]= count
  clusters[i,3]= length(clustspecies)
  i = i+1
}

orfancount = 0
orfan = c()
singcopycount = 0
singcopy = c()
multicopycount = 0
multicopy = c()
for (i in seq(1,nrow(clusters)))
{
  if (as.numeric(clusters[i,3]) <= 0.3 * numspecies)
  {
    orfancount = orfancount + 1
    orfan = c(orfan, as.numeric(clusters[i,2]))
  }
  else if (as.numeric(clusters[i,2]) <= 1.3 * numspecies)
  {
    singcopycount = singcopycount + 1
    singcopy = c(singcopy, as.numeric(clusters[i,2]))
  }
  else
  {
    multicopycount = multicopycount + 1
    multicopy = c(multicopy, as.numeric(clusters[i,2]))
  }
}
result = c(orfancount, singcopycount, multicopycount)

pdf("../results/cluster-size-dist.pdf")
m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
temp <- split.screen(m)
mar.default <- c(5,4,4,2) + 0.1
cols = c("green1","orange","gray35")
screen(1)
par(mar = mar.default + c(0, 1, 0, 0))
barplot(result, col=cols, yaxt="n", ylab='Frequency', cex.lab=1.5)
axis(1, at=c(0.7,1.9,3.1), labels=c('ORFan', 'Single-copy', 'Multiple-copy'), cex.axis=1.5)
axis(2, at=seq(0,4000,1000), labels=c(0,NA,NA,NA,4000), cex.axis=1.5)
screen(2)
par(mar = mar.default + c(0, 1, 0, 0))
barplot(table(factor(orfan, levels=1:100)), col=cols[1], ylim=c(1,53), horiz=TRUE, yaxt="n", xpd=FALSE, xaxt="n", xlab='Frequency',
        cex.lab=1.5, ylab='Cluster size')
axis(1, at=seq(0,2000,500), labels=c(0,NA,NA,NA,2000), cex.axis=1.5)
axis(2, at=c(0,13,26,39,52), labels=c(0,NA,26,NA,52), cex.axis=1.5)
screen(3)
barplot(table(factor(singcopy, levels=1:100)), col=cols[2], ylim=c(1,52), horiz=TRUE, yaxt="n", xpd=FALSE, xaxt="n", xlab='Frequency',
        cex.lab=1.5)
axis(1, at=seq(0,1000,200), labels=c(0,NA,NA,NA,NA,1000), cex.axis=1.5)
axis(2, at=c(0,13,26,39,52), labels=c(0,NA,26,NA,52), cex.axis=1.5)
screen(4)
barplot(table(factor(multicopy, levels=1:100)), col=cols[3], ylim=c(1,52), horiz=TRUE, yaxt="n", xpd=FALSE, xaxt="n", xlab='Frequency',
        cex.lab=1.5)
axis(1, at=seq(0,70,10), labels=c(0,NA,NA,NA,NA,NA,NA,70), cex.axis=1.5)
axis(2, at=c(0,13,26,39,52), labels=c(0,NA,26,NA,52), cex.axis=1.5)
close.screen(all.screens = TRUE)

dev.off()