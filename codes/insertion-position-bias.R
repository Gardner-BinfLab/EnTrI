insertion_positions <- read.table("../results/insertion-position-bias.out", row.names=1)
pdf("../results/insertion-position-bias.pdf")
toplot <- colMeans(insertion_positions[,1:ncol(insertion_positions)-1])
cols = c(rep("khaki4",60), rep("slategray",60), rep("palevioletred4",60))
midpoints <- barplot(toplot, main="All", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,180)], labels=c('5\'','3\''), cex.axis=1.5)
legend(61,0.7, c(expression(1 %->% 60),"internal", expression((L-60+1) %->% L)), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[61],cols[121]), bg="white")
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
axis(1, at=midpoints[c(1,180)], labels=c('5\'','3\''), cex.axis=1.5)
legend(61,0.5, c(expression(1 %->% 60),"internal", expression((L-60+1) %->% L)), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[61],cols[121]), bg="white")
toplot <- colMeans(nonessential)
midpoints <- barplot(toplot, main="Non essential", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,180)], labels=c('5\'','3\''), cex.axis=1.5)
legend(61,0.6, c(expression(1 %->% 60),"internal", expression((L-60+1) %->% L)), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[61],cols[121]), bg="white")
toplot <- colMeans(beneficialloss)
midpoints <- barplot(toplot, main="Beneficial loss", xaxt="n", xlab="position", ylab="mean ii", cex.axis=1.5, cex.lab=1.5, cex.main=2, col=cols, border=NA)
axis(1, at=midpoints[c(1,180)], labels=c('5\'','3\''), cex.axis=1.5)
legend(61,1.5, c(expression(1 %->% 60),"internal", expression((L-60+1) %->% L)), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[61],cols[121]), bg="white")
dev.off()