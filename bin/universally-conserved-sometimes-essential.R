library(gplots)
ucse <- read.csv("~/EnTrI/results/interesting_genes/universally-conserved_sometimes-essential.tsv", header=TRUE, sep='\t')
to.plot <- as.matrix(ucse[,3:18])
rownames(to.plot) = ucse[,1]
breaks <- seq(min(to.plot), max(to.plot), length.out=51)
my_palette <- colorRampPalette(c("cyan", "palevioletred4"))(n = 50)
pdf("~/EnTrI/figures/heatmap-ucse.pdf", width = 15, height = 30)
heatmap.2(to.plot, margins = c(10, 10), main = "Universally conserved sometimes essential", density.info="none", trace="none", keysize=1,
          lmat=rbind( c(0, 3), c(2,1), c(0,4)), lhei=c(0.03, 1, 0.05 ), lwid = c(0.05,1), Rowv=FALSE, Colv=FALSE, dendrogram='none',
          col=my_palette, breaks = breaks)
dev.off()