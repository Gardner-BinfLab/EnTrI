library(RColorBrewer)

dbpath <- '../results/core-ancestral/ancestral-not-core.heatmap'
db <- read.table(dbpath, stringsAsFactors = FALSE, sep='\t', header = TRUE)
sorted <- db[ order(db$path_name, db$Gene), ]
toplot <- sorted[,seq(3,16)]
toplot[toplot == 'X'] <- NA
toplot <- data.matrix(toplot)
rownames(toplot) = make.names(sorted$Gene, unique=TRUE)
annotation <- data.frame(row.names = make.names(sorted$Gene, unique=TRUE), category = sorted$path_name)

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

nHalf <- 50
Min <- min(toplot, na.rm=TRUE)
Max <- max(toplot, na.rm=TRUE)
Thresh <- 0
eps <- 1e-10
rc1 <- colorRampPalette(colors = c("palevioletred4", "pink"), space="Lab")(nHalf)
rc2 <- colorRampPalette(colors = c("paleturquoise", "lightsteelblue4"), space="Lab")(nHalf)
ramps <- c(rc1,rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf)
rb2 <- seq(Thresh+eps, Max, length.out=nHalf)
rampbreaks <- c(rb1, rb2)
col_fun = circlize::colorRamp2(rampbreaks, ramps)
# ha = HeatmapAnnotation(levels = rampbreaks[c(1,25,50,75,100)], col = list(NPEQ = ramps[c(1,25,50,75,100)]))

dups = sorted[,31:44]
rownames(dups) <- rownames(toplot)
colnames(dups) <- colnames(toplot)

inparalogs = sorted[,45:60]
rownames(inparalogs) <- rownames(toplot)
colnames(inparalogs) <- colnames(toplot)

pdf("~/EnTrI/figures/ancestral-not-core_heatmap.pdf", height = 60)
Heatmap(toplot, col = col_fun, name = 'log(ii/essentiality_threshold)', cluster_rows = FALSE, cluster_columns = FALSE, split = annotation, gap = unit(5, "mm"),
        na_col = 'black', heatmap_legend_param = list(color_bar = "continuous"),
        # combined_name_fun = function(x) 
        # {
        #   grid.text(x,rot=90)
        # },
        cell_fun = function(j, i, x, y, width, height, fill) 
        {
          if(dups[i, j] > 1 & !is.na(toplot[i, j])) 
          {
            grid.text("D", x = x, y = y)
          }
          # if(inparalogs[i, j] > 1 & !is.na(toplot[i, j])) 
          # {
          #   grid.text("IP", x = x, y = y)
          # }
          # else if(dups[i, j] > 1 & !is.na(toplot[i, j]))
          # {
          #   grid.text("OP", x = x, y = y)
          # }
        } # https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.htm
)
dev.off()