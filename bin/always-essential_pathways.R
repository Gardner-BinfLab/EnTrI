library(stringr)
dbpath <- "~/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat"
datapath <- "~/EnTrI/results/interesting_genes/universally-conserved_always-essential.tsv"
writepath <- "~/EnTrI/results/interesting_genes/universally-conserved_always-essential-with-modules.tsv"
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
all = length(unique(sorted$Gene))
pathways.names <- unique(sorted$path_name)
pathways <- rep(0, length(pathways.names))
names(pathways) <- pathways.names
for (i in seq(1, nrow(sorted)))
{
  len = sum(sorted$Gene==sorted$Gene[i])
  pathways[sorted$path_name[i]] = pathways[sorted$path_name[i]] + 1/len
}
percentage = (pathways/all)*100
percentage <- cbind(pathways.names, percentage)
row.names(percentage) <- NULL
colnames(percentage) <- c('pathway', 'percent')
percentage <- as.data.frame(percentage)
percentage$percent <- as.numeric(levels(percentage$percent))[percentage$percent]
# percentage <- percentage[ order(percentage$percent), ]
library(ggplot2)
pdf('~/EnTrI/figures/pathways-alwaysessential.pdf')
ggplot(data=percentage, aes(x= reorder(pathway, percent), y = percent)) + geom_bar(stat = 'identity') + coord_flip() + xlab("Pathways") +
  ylab("Percentage") + labs(title="Pathways in always essential genes")
dev.off()
