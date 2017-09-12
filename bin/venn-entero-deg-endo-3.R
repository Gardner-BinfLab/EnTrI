library('UpSetR')
subsets <- read.table('../results/venn-entero-deg-endo/venn.txt', sep = '\t')
expressionInput <- setNames(as.character(subsets$V2), subsets$V1)
pdf('../figures/venn-entero-deg-endo.pdf')
upset(fromExpression(expressionInput), order.by = "freq",
      sets=rev(c('Enterobacteriaceae', 'Endosymbionts', 'Gammaproteobacteria', 'Proteobacteria', 'All')), keep.order = TRUE,
      text.scale = 2)
dev.off()