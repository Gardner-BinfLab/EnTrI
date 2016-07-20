library(stringr)
clusters_path <- "../results/hieranoid/hieranoid-result.txt"
etec_path <- "../results/insertion-indices/normalised-insertion-indices/ETEC.txt"
k12_path <- "../results/ecogene-k12.txt"
eteciis <- as.matrix(read.table(etec_path))
rownames(eteciis) <- eteciis[,1]
eteciis <- eteciis[,-1]
# hist(as.numeric(iis[,2]), breaks=300)
# list_of_essential_genes["ETEC"] <- list(iis[,1][iis[,3]=='essential'])
biis <- as.matrix(read.table(k12_path))
cluster <- readLines(clusters_path)
numclusters = length(cluster)
k12conserved = c()
k12unconserved = c()
k12essential = c()
for (clusterindex in (1:numclusters))
{
  line = cluster[clusterindex]
  bmatches = str_match_all(line, 'b[[:digit:]]+')
  etecmatches = str_match_all(line, 'ETEC_[[:digit:]]+')
  if (is.na(bmatches[[1]][1]) & !is.na(etecmatches[[1]][1]))
  {
    k12unconserved = c(k12unconserved, etecmatches[[1]])
  }
  else if (!is.na(etecmatches[[1]][1]))
  {
    essentiality = 0
    for (i in (1:nrow(bmatches[[1]])))
    {
      if (bmatches[[1]][i] %in% biis)
      {
        essentiality = 1
      }
      if (essentiality)
      {
        k12essential = c(k12essential, etecmatches[[1]])
      }
      else
      {
        k12conserved = c(k12conserved, etecmatches[[1]])
      }
    }
  }
}
for (i in (length(k12conserved):1))
{
  if (!(k12conserved[i] %in% row.names(eteciis)))
  {
    k12conserved <- k12conserved[-i]
  }
}
etecessential_k12unconserved = c()
for (i in (length(k12unconserved):1))
{
  if (!(k12unconserved[i] %in% row.names(eteciis)))
  {
    k12unconserved <- k12unconserved[-i]
  }
  else if(eteciis[k12unconserved[i],1] <= 0.2799)
  {
    etecessential_k12unconserved = c(etecessential_k12unconserved, k12unconserved[i])
  }
}
etecinessential_k12essential = c()
for (i in (length(k12essential):1))
{
  if (!(k12essential[i] %in% row.names(eteciis)))
  {
    k12essential <- k12essential[-i]
  }
  else if(eteciis[k12essential[i],1] > 0.1324)
  {
    etecinessential_k12essential = c(etecinessential_k12essential, k12essential[i])
  }
}
# pdf("../figures/ETEC-vs-K12.pdf")
h1=hist(as.numeric(eteciis[k12essential,1]), breaks=seq(0,5,0.03), plot=F)$counts
h2=hist(as.numeric(eteciis[k12conserved,1]), breaks=seq(0,5,0.03), plot=F)$counts
h3=hist(as.numeric(eteciis[k12unconserved,1]), breaks=seq(0,5,0.03), plot=F)$counts
barplot(rbind(h1,h2,h3),col=c("darkgoldenrod4","turquoise4", "darkmagenta"), main="Escherichia coli ETEC ETEC")
legend(50,180, c("essential in k-12","conserved non-essential in k-12", "not conserved in k-12"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
       col=c("darkgoldenrod4","turquoise4", "darkmagenta"), bty="n")
lines(c(0.1324/0.03, 0.1324/0.03), c(0,300), col="red")
# dev.off()
