library(ComplexHeatmap)
library(stringr)
toplot <- read.table('~/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap.tsv', header = TRUE)
rownames(toplot) = make.names(toplot$gene, unique=TRUE)
toplot <- toplot[,seq(2,6)]
toplot <- toplot[ order(toplot$Citrobacter, toplot$Salmonella, toplot$Escherichia, toplot$Klebsiella, toplot$Enterobacter,
                        row.names(toplot), decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), method = 'radix'), ]
pdf("~/EnTrI/figures/accessory_essential_enterobacteriaceae_heatmap.pdf", height = 75, width = 4)
Heatmap(toplot, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = 'gray'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()