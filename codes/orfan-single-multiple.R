library(stringr)
list_of_files <- list.files(path="../results/homclust/EFam-clusters", pattern="*.txt", full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "b")
genuses = c(1, 2, 3, 2, 2, 3, 3, 4, 4, 4, 4, 4, 2)
names(genuses) <- names
numspecies = length(names)
clusters = matrix(nrow=length(list_of_files), ncol=5)
i = 1
for (filename in list_of_files)
{
  clust_tbl =  read.table(filename, header = FALSE)
  clustspecies = c()
  for (item in clust_tbl[,2])
  {
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
  unique_clustspecies = unique(clustspecies)
  clusters[i,1]= basename(filename)
  clusters[i,2]= nrow(clust_tbl)
  clusters[i,3]= length(unique_clustspecies)
  clusters[i,4]= length(table(clustspecies)[table(clustspecies)>1])
  cluster_genuses = c()
  for (item in unique_clustspecies)
  {
    cluster_genuses <- c(cluster_genuses, genuses[unique_clustspecies])
  }
  clusters[i,5]= length(unique(cluster_genuses))
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
  if (as.numeric(clusters[i,5]) <= 1)
  {
    orfancount = orfancount + 1
    orfan = c(orfan, as.numeric(clusters[i,2]))
  }
  else if (as.numeric(clusters[i,3]) - as.numeric(clusters[i,4]) >= 0.7 * as.numeric(clusters[i,3]))
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

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0)) 
n= 30
cols = c("green1","orange", "gray35")
midpoints <- barplot(rbind(table(factor(orfan, levels=1:100)),table(factor(singcopy, levels=1:100)),
                           table(factor(multicopy, levels=1:100))), col = cols, xlim=c(0,n), xaxt="n", border = NA,
                     xlab="Cluster size",ylab="Frequency", main ="Cluster size distribution", , cex.lab = 2, cex.axis = 1.5, cex.main = 2)
axis(1, at=midpoints[seq(2,n,2)], labels=c(rbind(seq(2,n,4),NA))[1:(n/2)], cex.axis=1.5)
legend(7,2300, c("Genus Specific","Single copy", "Multiple copy"), lty=c(1,1,1),cex=1.5, bty="n", lwd=c(4,4, 4),col=cols)

dev.off()