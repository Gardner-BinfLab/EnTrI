library(pheatmap)
library(RColorBrewer)
library(stringr)
dbpath <- "~/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat"
datapath <- "~/EnTrI/results/interesting_genes/sometimes-essential-marked-dup.tsv"
writepath <- "~/EnTrI/results/interesting_genes/sometimes-essential-marked-dup-with-modules.tsv"
db <- read.csv(dbpath, sep='\t', stringsAsFactors = FALSE)
for (i in seq(1,nrow(db)))
{
  db$path_name[i] = str_match(db$path_name[i], '([[:print:]]+)- Escherichia coli K\\-12 MG1655')[2]
  
}
data <- read.csv(datapath, sep='\t', stringsAsFactors = FALSE)
merged <- merge(x=data, y=db, by.x="EsCoMG1655.locus.tag", by.y="gene_id")
merged <- merged[,c(seq(2,34),1,35,36)]
write.table(merged, file=writepath, quote=FALSE, sep = '\t', row.names = FALSE)
sorted <- merged[ order(merged$path_name, merged$Gene), ]
toplot <- sorted[,seq(3,18)]
toplot$EsCoMG1655.NPEQ[toplot$EsCoMG1655.NPEQ == 7] <- 9
toplot[toplot == 'X'] <- NA
toplot <- data.matrix(toplot)
rownames(toplot) = make.names(sorted$Gene, unique=TRUE)
annotation <- data.frame(row.names = make.names(sorted$Gene, unique=TRUE), category = sorted$path_name)
gaps = c()
for (i in seq(2, length(annotation$category)))
  if (annotation$category[i] != annotation$category[i-1])
    gaps = c(gaps, i-1)
my_palette <- c(colorRampPalette(c("cyan","gray"))(n = 15),colorRampPalette(c("gray","palevioletred4"))(n = 35))
newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotation$category))))
mycolors <- newCols(length(unique(annotation$category)))
names(mycolors) <- unique(annotation$category)
mycolors <- list(category = mycolors)

pheatmap(toplot, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = annotation, gaps_row = gaps, cellheight = 10, cellwidth = 20,
         file = "~/EnTrI/figures/interesting-genes_heatmap.pdf", color = my_palette, annotation_colors = mycolors)


# library(ggplot2)
# library(stringr)
# dbs <- "~/EnTrI/results/KEGG_modules/"
# data <- "~/EnTrI/results/interesting_genes/sometimes-essential-no-dup.tsv"
# list_of_files <- list.files(path=dbs, full.names=T, recursive=FALSE)
# db=matrix(,nrow=0,ncol=3)
# for (filename in list_of_files)
# {
#   dbfile = as.matrix(read.csv(filename, as.is=TRUE, header=TRUE, sep="\t"))
#   # for (i in seq(1,nrow(dbfile)))
#   # {
#   #   dbfile[i,3] = str_match(dbfile[i,3], '([[:print:]]+)-[[ Citrobacter rodentium]]|[[ Enterobacter cloacae subsp. cloacae NCTC 9394]]|
#   #                           [[ Escherichia coli O78:H11:K80 H10407 (ETEC)]]|[[ Escherichia coli K\\-12 MG1655]]|
#   #                           [[ Salmonella enterica subsp. enterica serovar Enteritidis P125109]]|
#   #                           [[ Salmonella enterica subsp. enterica serovar Typhimurium D23580]]|
#   #                           [[ Salmonella enterica subsp. enterica serovar Typhimurium SL1344]]|
#   #                           [[ Salmonella enterica subsp. enterica serovar Typhi Ty2]]|
#   #                           [[ Escherichia coli O25b:K100:H4-ST131 EC958 (UPEC)]]')[2]
#   #   
#   # }
#   db = rbind(db,dbfile)
# }
# sometimesess <- read.csv(data, header = TRUE, sep='\t')
# sometimesess <- sometimesess[,c("CiRoICC168.locus.tag", "EnClNCTC9394.locus.tag", "EsCoH10407.locus.tag", "EsCoMG1655.locus.tag",
#                                 "SaEnP125109.locus.tag", "SaTyD23580.locus.tag", "SaTySL1344.locus.tag", "SaTyTy2.locus.tag", "EsCoEC958.locus.tag")]
# genes <- unlist(sometimesess)
# genes <- genes[genes != "X"]
# 
# enrichment = data.frame(matrix(nrow=length(unique(db[,3])),ncol=4),row.names =unique(db[,3]) )
# for (i in seq(1,nrow(enrichment)))
# {
#   for (j in seq(1,ncol(enrichment)))
#   {
#     enrichment[i,j]=0
#   }
# }
# 
# len_genes = 0
# for (item in genes)
# {
#   for (i in seq(1,nrow(db)))
#   {
#     if (item == db[i,1])
#     {
#       enrichment[db[i,3], 1] = enrichment[db[i,3], 1] + 1
#       len_genes = len_genes + 1
#     }
#   }
# }
# 
# for (i in seq(1,nrow(db)))
# {
#   enrichment[db[i,3], 2] = enrichment[db[i,3], 2] + 1
# }
# 
# len_db = length(db[,1])
# for (i in seq(1,nrow(enrichment)))
# {
#   # fisherresult = fisher.test(matrix(c(enrichment[i,1],enrichment[i,2],(len_genes-enrichment[i,1]), (len_db-enrichment[i,2])),nrow=2,byrow=TRUE), alternative = "greater")
#   # enrichment[i,3] = fisherresult$p.value
#   
#   x = enrichment[i,1]
#   m = enrichment[i,2]
#   n = len_db - enrichment[i,2]
#   k = len_genes
#   enrichment[i,3] = phyper(x,m,n,k, lower.tail = FALSE)
#   
#   
# }
# enrichment[,4] = p.adjust(enrichment[,3], method = "BH")
# enrichment <- enrichment[order(enrichment[,4], enrichment[,3]),]
# 
# pathways=data.frame()
# i=1
# while(enrichment[i,4] < 0.05)
# {
#   j = 1
#   while(db[j,3] != row.names(enrichment)[i])
#   {
#     j = j + 1
#   }
#   pathways[i,1] = db[j,3]
#   pathways[i,2] = -log10(enrichment[i,4])
#   i = i + 1
# }
# colnames(pathways) = c('word', 'pvalue')
# pathways = pathways[1:20,]
# 
# pdf("../figures/sometimes-essential-modules.pdf")
# ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
#   geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Sometimes essential modules") +
#   theme(text = element_text(size=15))
# dev.off()