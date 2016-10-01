library(ggplot2)
db <- as.matrix(read.table("~/EnTrI/results/KEGG/cirobacter_rodentium_ICC168.dat", header=TRUE))
genes <- as.matrix(read.table("~/EnTrI/results/KEGG/ROD-essential.txt", header=TRUE))
enrichment = data.frame(matrix(nrow=length(unique(db[,2])),ncol=4),row.names =unique(db[,2]) )
for (i in seq(1,nrow(enrichment)))
{
  for (j in seq(1,ncol(enrichment)))
  {
    enrichment[i,j]=0
  }
}
len_genes = 0
for (item in genes)
{
  for (i in seq(1,nrow(db)))
  {
    if (item == db[i,1])
    {
      enrichment[db[i,2], 1] = enrichment[db[i,2], 1] + 1
      len_genes = len_genes + 1
    }
  }
}

for (i in seq(1,nrow(db)))
{
  enrichment[db[i,2], 2] = enrichment[db[i,2], 2] + 1
}

len_db = length(db[,1])
for (i in seq(1,nrow(enrichment)))
{
  fisherresult = fisher.test(matrix(c(enrichment[i,1],enrichment[i,2],(len_genes-enrichment[i,1]), (len_db-enrichment[i,2])),nrow=2,byrow=TRUE), alternative = "greater")
  enrichment[i,3] = fisherresult$p.value
}
enrichment[,4] = p.adjust(enrichment[,3], method = "BH")
enrichment <- enrichment[order(enrichment[,4], enrichment[,3]),]

pathways=data.frame()
i=1
while(enrichment[i,4] < 0.05)
{
  j = 1
  while(db[j,2] != row.names(enrichment)[i])
  {
    j = j + 1
  }
  pathwayname = strsplit(db[j,3],"-")[[1]][1]
  pathways[i,1] = pathwayname
  pathways[i,2] = -log10(enrichment[i,4])
  i = i + 1
}
colnames(pathways) = c('word', 'pvalue')

# pdf("../figures/pathways.pdf")
ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Essential genes") +
  theme(text = element_text(size=15))
# dev.off()