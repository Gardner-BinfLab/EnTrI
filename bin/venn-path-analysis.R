library(ggplot2)
library(stringr)
dbpath <- "~/EnTrI/results/homologs.path"
essentiality <- "~/EnTrI/results/venn-entero-deg-endo-del/entero1deg0endo1.txt"
db = as.matrix(read.table(dbpath, as.is=TRUE, header=TRUE, sep="\t"))
for (i in seq(1,nrow(db)))
{
  db[i,1] = paste('entero', db[i,1], sep = '-')
  db[i,3] = str_match(db[i,3], '([[:print:]]+) - Escherichia coli K\\-12 MG1655')[2]
}
essentialityfile = as.matrix(read.table(essentiality, as.is=TRUE))
genes = c()
firstrow = essentialityfile[1,]
num = which(startsWith(firstrow,'entero'))
for (i in seq(1,nrow(essentialityfile)))
{
  genes=c(genes,essentialityfile[i,num])
}

enrichment = data.frame(matrix(nrow=length(unique(db[,3])),ncol=4),row.names =unique(db[,3]) )
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
      enrichment[db[i,3], 1] = enrichment[db[i,3], 1] + 1
      len_genes = len_genes + 1
    }
  }
}

for (i in seq(1,nrow(db)))
{
  enrichment[db[i,3], 2] = enrichment[db[i,3], 2] + 1
}

len_db = length(db[,1])
for (i in seq(1,nrow(enrichment)))
{
  # fisherresult = fisher.test(matrix(c(enrichment[i,1],enrichment[i,2],(len_genes-enrichment[i,1]), (len_db-enrichment[i,2])),nrow=2,byrow=TRUE), alternative = "greater")
  # enrichment[i,3] = fisherresult$p.value
  
  x = enrichment[i,1]
  m = enrichment[i,2]
  n = len_db - enrichment[i,2]
  k = len_genes
  enrichment[i,3] = phyper(x,m,n,k, lower.tail = FALSE)
}
enrichment[,4] = p.adjust(enrichment[,3], method = "BH")
enrichment <- enrichment[order(enrichment[,4], enrichment[,3]),]

pathways=data.frame()
i=1
while(enrichment[i,4] < 0.05)
{
  j = 1
  while(db[j,3] != row.names(enrichment)[i])
  {
    j = j + 1
  }
  pathways[i,1] = db[j,3]
  pathways[i,2] = -log10(enrichment[i,4])
  i = i + 1
}
colnames(pathways) = c('word', 'pvalue')
pathways = pathways[1:20,]
pathways = na.omit(pathways)
pdf("../figures/entero1deg0endo1.pdf")
ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("(Entero & Endo) - DEG") +
  theme(text = element_text(size=15))
dev.off()
