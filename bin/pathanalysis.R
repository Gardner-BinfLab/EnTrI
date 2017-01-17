library(ggplot2)
library(stringr)
dbs <- "~/EnTrI/results/KEGG/"
essentiality <- "~/EnTrI/results/biases/normalised-pca/"
list_of_files <- list.files(path=dbs, full.names=T, recursive=FALSE)
db=matrix(,nrow=0,ncol=3)
for (filename in list_of_files)
{
  dbfile = as.matrix(read.table(filename, as.is=TRUE, header=TRUE, sep="\t"))
  for (i in seq(1,nrow(dbfile)))
  {
    dbfile[i,3] = str_match(dbfile[i,3], '([[:print:]]+)-[[ Citrobacter rodentium]]|[[ Enterobacter cloacae subsp. cloacae NCTC 9394]]|
                            [[ Escherichia coli O78:H11:K80 H10407 (ETEC)]]|[[ Escherichia coli K\\-12 MG1655]]|
                            [[ Salmonella enterica subsp. enterica serovar Enteritidis P125109]]|
                            [[ Salmonella enterica subsp. enterica serovar Typhimurium D23580]]|
                            [[ Salmonella enterica subsp. enterica serovar Typhimurium SL1344]]|
                            [[ Salmonella enterica subsp. enterica serovar Typhi Ty2]]|
                            [[ Escherichia coli O25b:K100:H4-ST131 EC958 (UPEC)]]')[2]
    
  }
  db = rbind(db,dbfile)
}
list_of_files <- list.files(path=essentiality, full.names=T, recursive=FALSE)
genes = c()
for (filename in list_of_files)
{
  essentialityfile = as.matrix(read.table(filename, as.is=TRUE, header=TRUE))
  for (i in seq(1,nrow(essentialityfile)))
  {
    if (essentialityfile[i,3] == "essential")
    {
      genes=c(genes,essentialityfile[i,1])
    }
  }
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

pdf("../figures/essential-pathways.pdf")
ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Essential genes") +
  theme(text = element_text(size=15))
dev.off()