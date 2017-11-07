library(ComplexHeatmap)
toplot <- read.table('~/EnTrI/results/venn-entero-deg-endo/heatmap.tsv', header = TRUE)
rownames(toplot) = make.names(toplot$gene, unique=TRUE)
toplot <- toplot[,seq(2,6)]
toplot <- toplot[ order(toplot$Enterobacteriaceae, toplot$Endosymbiont, toplot$Gammaproteobacteria, toplot$Proteobacteria, toplot$Bacteria,
                        row.names(toplot), decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), method = 'radix'), ]
pdf("~/EnTrI/figures/essential_genes_heatmap.pdf", height = 60, width = 3)
Heatmap(toplot, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = 'black'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()
