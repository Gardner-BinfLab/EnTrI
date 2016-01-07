biasespath <- "../results/check-biases.out"
biasestable = read.table(biasespath, header = FALSE)
colnames(biasestable) <- c("name", "ii", "dist", "gc")

pdf("../results/biases.pdf")

plot(biasestable$dist, biasestable$ii, pch = '.', ylim=c(0,3))
lines(lowess(biasestable$dist, biasestable$ii, f = 0.2), col = 2, lwd=5)
regln <- lm(biasestable$dist ~ biasestable$ii)
abline(regln, col=3, lwd=5)

plot(biasestable$gc, biasestable$ii, pch = '.', ylim=c(0,3))
lines(lowess(biasestable$gc, biasestable$ii, f = 0.2), col = 2, lwd=5)
regln <- lm(biasestable$gc ~ biasestable$ii)
abline(regln, col=3, lwd=5)

dev.off()