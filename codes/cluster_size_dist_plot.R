list_of_files <- list.files(path="../results/homclust/EFam-clusters", pattern="*.txt", full.names=T, recursive=FALSE)
count = c()
clusters = matrix(nrow=length(list_of_files), ncol=2)
i = 1
for (filename in list_of_files)
{
  count = c(count, length(readLines(filename)))
  clusters[i,1]= basename(filename)
  clusters[i,2]=count[i]
  i = i+1
}
count_table = table(count)

pdf("../results/cluster-size-dist.pdf")

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0)) 
n= 30
cols = c(rep("green1", 5), rep("orange", 13), rep("gray35", 96))
midpoints <- barplot(count_table,xlab="Cluster size",ylab="Frequency", main ="Cluster size distribution", cex.lab = 2, cex.axis = 1.5, cex.main = 2,
        #col = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1),
        col=cols,
        xlim=c(0,30), xaxt="n"
        #, yaxt="n"
        )
#axis(2, at=c(0,10,20,30,40), labels=c("0","100","400","900","1600"))
axis(1, at=midpoints[seq(2,n,2)], labels=seq(2,n,2))
#plot(ecdf(count), xlim=c(0,20))
lines(c(6.1, 6.1), c(-100,3000), col = "red", lwd=3, lty = 2)
lines(c(21.7, 21.7), c(-100,3000), col = "red", lwd=3, lty = 2)
legend(7,2300, c("ORFan","Single copy", "Multiple copy"), lty=c(1,1,1),cex=1.5, bty="n", lwd=c(4,4, 4),col=c("green1","orange", "gray35"))

dev.off()