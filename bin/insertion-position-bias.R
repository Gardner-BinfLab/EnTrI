insertion_positions <- read.table("../results/insertion-indices/insertion-position-bias.out", row.names=1)
pdf("../figures/insertion-position-bias.pdf")
toplot <- colMeans(insertion_positions[,1:ncol(insertion_positions)-1])
cols = c(rep("khaki4",5), rep("slategray",75), rep("palevioletred4",20))
midpoints <- barplot(toplot, main="All", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=1.5)
legend(41,0.7, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
essential = c()
nonessential = c()
beneficialloss = c()
for (i in seq(1,nrow(insertion_positions)))
{
  if (insertion_positions[i,ncol(insertion_positions)] < 0.2)
  {
    essential = rbind(essential, insertion_positions[i,1:ncol(insertion_positions)-1]) 
  }
  else if (insertion_positions[i,ncol(insertion_positions)] < 2)
  {
    nonessential = rbind(nonessential, insertion_positions[i,1:ncol(insertion_positions)-1])
  }
  else
  {
    beneficialloss = rbind(beneficialloss, insertion_positions[i,1:ncol(insertion_positions)-1]) 
  }
}
toplot <- colMeans(essential)
midpoints <- barplot(toplot, main="Essential", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=1.5)
legend(41,0.3, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
toplot <- colMeans(nonessential)
midpoints <- barplot(toplot, main="Non essential", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=1.5)
legend(41,0.6, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
toplot <- colMeans(beneficialloss)
midpoints <- barplot(toplot, main="Beneficial loss", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=1.5)
legend(41,1.5, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
dev.off()