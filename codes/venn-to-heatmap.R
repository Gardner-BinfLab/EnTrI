library("gplots")
venspath <- "../results/pairwise-venn-diagrams/"
list_of_files <- list.files(path=venspath, full.names=T, recursive=FALSE)
plot_mains = c("Ess(rows)-Abs(columns)", "Ess-Ess", "Ess(rows)-NoEs(columns)", "Pres(rows)-Abs(columns)", "Pres-Pres")
i = 1

my_palette <- colorRampPalette(c("steelblue4", "slategray2"))(n = 100)

pdf("../results/essentiality-heatmap.pdf")

for (filename in list_of_files)
{
  heattable <- as.matrix(read.table(filename))
  names <- heattable[1,2:nrow(heattable)]
  heattable <- heattable[-1,]
  heattable <- heattable[,-1]
  rownames(heattable) <- names
  colnames(heattable) <- names
  heattable <- heattable[c(5,6,7,8,9,2,3,4,1,10,11,12),c(5,6,7,8,9,2,3,4,1,10,11,12)]
  class(heattable) <- "numeric" 
  heatmap.2(heattable, margins = c(15, 15), main = plot_mains[i], density.info="none", trace="none", keysize=1,
            lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 6, 2 ), Rowv=FALSE, Colv=FALSE, dendrogram='none', col=my_palette)
  i = i+1
}

dev.off()