library("gplots")
venspath <- "../results/pairwise-venn-diagrams/"
list_of_files <- list.files(path=venspath, full.names=T, recursive=FALSE)
plot_mains = c("Rows: essential, Columns: absent", "Rows: essential, Columns: essential", "Rows: essential, Columns: non-essential",
               "Rows: present, Columns: absent", "Rows: present, Columns: present")
i = 1

pdf("../results/essentiality-heatmap.pdf")

for (filename in list_of_files)
{
  heattable <- as.matrix(read.table(filename))
  names <- heattable[1,2:nrow(heattable)]
  heattable <- heattable[-1,]
  heattable <- heattable[,-1]
  rownames(heattable) <- names
  colnames(heattable) <- names
  class(heattable) <- "numeric" 
  heatmap.2(heattable, margins = c(15, 15), main = plot_mains[i], density.info="none", trace="none", keysize=1,
            lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 6, 2 ))
  i = i+1
}

dev.off()