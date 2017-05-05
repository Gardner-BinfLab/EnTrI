library("dbscan")
library(ComplexHeatmap)
datapath <- "~/EnTrI/results/interesting_genes/universally-conserved.tsv"
data <- read.table(datapath, header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", comment.char = "")
data <- data[-c(143,1316,1360),]
names(data) <- gsub("NPEQ", "II", names(data))
iis <- data.matrix(data[,3:18])

data2path <- "~/EnTrI/results/interesting_genes/universally-conserved-log-ii-thresh.tsv"
data2 <- read.table(data2path, header = TRUE, sep = '\t', stringsAsFactors = FALSE, quote = "", comment.char = "")
data2 <- data2[-c(143,1316,1360),]
names(data2) <- gsub("NPEQ", "II", names(data2))
logiis <- data.matrix(data2[,3:18])

res <- dbscan(iis, minPts = 100, eps = 0.2)
# pdf('~/EnTrI/figures/conserved-genes-clustered.pdf')
plot(iis[,9],iis[,15],col=res$cluster+1, pch=20, xlab = names(data)[11], ylab = names(data)[17])
# dev.off()

# kmeansres <- kmeans(iis,2)
# plot(iis, col=kmeansres$cluster, pch=20, xlab = names(data)[11], ylab = names(data)[17])
# 
hieres <- hclust(dist(iis))
clusterCut <- cutree(hieres, 2)
# pdf('~/EnTrI/figures/conserved-genes-clustered.pdf')
plot(iis, col=clusterCut, pch=20)
# dev.off()

datatoplot <- cbind(data2[,c(1,seq(3,18))], res$cluster)
datatoplot <- datatoplot[ order(datatoplot$`res$cluster`, datatoplot$Gene), ]
annotation <- data.frame(row.names = make.names(datatoplot$Gene, unique=TRUE), category = c('Not always essential', 'Core essential')[datatoplot$`res$cluster`+1])
toplot <- data.matrix(datatoplot[,seq(2,17)])
rownames(toplot) <- datatoplot$Gene

nHalf <- 50
Min <- min(logiis)
Max <- max(logiis)
Thresh <- 0
eps <- 1e-10
rc1 <- colorRampPalette(colors = c("palevioletred4", "pink"), space="Lab")(nHalf)
rc2 <- colorRampPalette(colors = c("paleturquoise", "lightsteelblue4"), space="Lab")(nHalf)
ramps <- c(rc1,rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf)
rb2 <- seq(Thresh+eps, Max, length.out=nHalf)
rampbreaks <- c(rb1, rb2)
col_fun = circlize::colorRamp2(rampbreaks, ramps)

pdf("~/EnTrI/figures/conserved-genes-clustered.pdf", height = 300)
Heatmap(toplot, col = col_fun, name = 'log(ii/essentiality_threshold)', cluster_rows = FALSE, cluster_columns = FALSE, split = annotation, gap = unit(5, "mm"),
        na_col = 'black', heatmap_legend_param = list(color_bar = "continuous"))
dev.off()