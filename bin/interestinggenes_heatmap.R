# library(pheatmap)
library(RColorBrewer)
library(stringr)
dbpath <- "~/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat"
datapath <- "~/EnTrI/results/interesting_genes/sometimes-essential-marked-dup.tsv"
writepath <- "~/EnTrI/results/interesting_genes/sometimes-essential-marked-dup-with-pathways.tsv"
# datapath <- "~/EnTrI/results/interesting_genes/universally-conserved_always-essential.tsv"
# writepath <- "~/EnTrI/results/interesting_genes/universally-conserved_always-essential-with-modules.tsv"
db <- read.csv(dbpath, sep='\t', stringsAsFactors = FALSE)
for (i in seq(1,nrow(db)))
{
  db$path_name[i] = str_match(db$path_name[i], '([[:print:]]+)- Escherichia coli K\\-12 MG1655')[2]

}
data <- read.csv(datapath, sep='\t', stringsAsFactors = FALSE, quote = "", comment.char = "")

#add always essential genes to the plot
datapath <- "~/EnTrI/results/interesting_genes/universally-conserved_always-essential-marked-dup.tsv"
data <- rbind(data, read.csv(datapath, sep='\t', stringsAsFactors = FALSE, quote = "", comment.char = ""))
names(data) <- gsub("NPEQ", "II", names(data))

merged <- merge(x=data, y=db, by.x="EsCoMG1655.locus.tag", by.y="gene_id")
# merged <- merged[,c(seq(2,34),1,35,36)]
merged <- merged[,c(seq(2,30),1,seq(31,60))]
write.table(merged, file=writepath, quote=FALSE, sep = '\t', row.names = FALSE)
sorted <- merged[ order(merged$path_name, merged$Gene), ]
toplot <- sorted[,seq(3,16)]
# toplot$EsCoMG1655.NPEQ[toplot$EsCoMG1655.NPEQ == 7] <- 9
toplot[toplot == 'X'] <- NA
toplot <- data.matrix(toplot)
rownames(toplot) = make.names(sorted$Gene, unique=TRUE)
annotation <- data.frame(row.names = make.names(sorted$Gene, unique=TRUE), category = sorted$path_name)
# gaps = c()
# for (i in seq(2, length(annotation$category)))
#   if (annotation$category[i] != annotation$category[i-1])
#     gaps = c(gaps, i-1)
# my_palette <- c(colorRampPalette(c("cyan","blue"))(n = 15),colorRampPalette(c("gray","palevioletred4"))(n = 35))
# newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation$category))))
# mycolors <- newCols(length(unique(annotation$category)))
# names(mycolors) <- unique(annotation$category)
# mycolors <- list(category = mycolors)
# pheatmap(toplot, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = annotation, gaps_row = gaps, cellheight = 10, cellwidth = 20,
#          file = "~/EnTrI/figures/interesting-genes_heatmap.pdf", color = my_palette, annotation_colors = mycolors)

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

toplot <- toplot[,c(4,10,9,5,6,7,8,12,11,13,14,3,2,1)]
colnames(toplot) <- c('Citrobacter.rodentium.ICC168', 'Salmonella.Typhi.Ty2', 'Salmonella.Enteritidis.P125109', 'Salmonella.Typhimurium.SL1344',
                      'Salmonella.Typhimurium.SL3261', 'Salmonella.Typhimurium.D23580', 'Salmonella.Typhimurium.A130',
                      'Escherichia.coli.UPEC.ST131', 'Escherichia.coli.ST131.EC958', 'Escherichia.coli.BW25113', 'Escherichia.coli.K-12.MG1655',
                      'Klebsiella.pneumoniae.RH201207', 'Klebsiella.pneumoniae.Ecl8', 'Enterobacter.cloacae.NCTC.9394')

# pdf("~/EnTrI/figures/interesting-genes-modules_heatmap.pdf", height = 70)
pdf("~/EnTrI/figures/interesting-genes_heatmap.pdf", height = 180)
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

ess = data[,3:16]
ess[ess == 'X'] <- NA
repetition = data[,31:44]

essdup <- sum(repetition > 1 & ess > Thresh, na.rm = TRUE)
nessdup <- sum(repetition > 1 & ess <= Thresh, na.rm = TRUE)
essndup <- sum(repetition == 1 & ess > Thresh, na.rm = TRUE)
nessndup <- sum(repetition == 1 & ess <= Thresh, na.rm = TRUE)
duplication.essentiality <- matrix(c(essdup,nessdup,essndup,nessndup), nrow = 2)
fisher.test(duplication.essentiality, alternative = "l")
