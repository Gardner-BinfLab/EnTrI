biasespath <- "../results/check-biases.out"
biasestable = read.table(biasespath, header = FALSE)
colnames(biasestable) <- c("name", "ii", "dist", "gc")

pdf("../results/biases.pdf")

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

plot(biasestable$dist, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "Distance from the origin", ylab = "insertion index",
     main = "Distance bias", cex.lab = 2, cex.axis = 2, cex.main =2)
lines(lowess(biasestable$dist, biasestable$ii, f = 0.2), col = 2, lwd=5)
legend(0.26,2.8,c("regression line", "lowess"), col = c(3,2), lty=1, lwd=5, cex = 1.5)
regln <- lm(biasestable$dist ~ biasestable$ii)
abline(regln, col=3, lwd=5)

plot(biasestable$gc, biasestable$ii, pch = '.', ylim=c(0,3), xlab = "GC content", ylab = "insertion index", main = "GC bias",
     cex.lab = 2, cex.axis = 2, cex.main =2)
# lines(lowess(biasestable$gc, biasestable$ii, f = 0.2), col = 2, lwd=5)
regln <- lm(biasestable$gc ~ biasestable$ii)
abline(regln, col=3, lwd=5)

dev.off()