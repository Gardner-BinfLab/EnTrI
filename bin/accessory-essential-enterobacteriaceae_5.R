library(ggplot2)
library(stringr)
dbpath <- "~/EnTrI/results/KEGG/escherichia_coli_K-12_MG1655.dat"
heatmappath <- "~/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap.tsv"
keggpath <- "~/EnTrI/results/define-core-accessory-hieranoid-fitch/heatmap-kegg.tsv"

db=matrix(,nrow=0,ncol=3)
dbfile = as.matrix(read.table(dbpath, as.is=TRUE, header=TRUE, sep="\t"))
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

heatmap = read.table(heatmappath, as.is=TRUE, header=TRUE, sep="\t")
kegg = read.table(keggpath, as.is=TRUE, header=FALSE, sep="\t")

ind = 2
start = 1
while(ind <= nrow(heatmap))
{
  row = heatmap[ind,2:6]
  prev_row = heatmap[ind-1,2:6]
  if(sum(row != prev_row))
  {
    genes = kegg[start:(ind-1),1]
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
        if (!is.na(item) & item == db[i,1])
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
      fisherresult = fisher.test(matrix(c(enrichment[i,1],enrichment[i,2],(len_genes-enrichment[i,1]), (len_db-enrichment[i,2])),nrow=2,byrow=TRUE), alternative = "greater")
      enrichment[i,3] = fisherresult$p.value
      
      # x = enrichment[i,1]
      # m = enrichment[i,2]
      # n = len_db - enrichment[i,2]
      # k = len_genes
      # enrichment[i,3] = phyper(x,m,n,k, lower.tail = TRUE)
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
    if (nrow(pathways) > 0)
    {
      colnames(pathways) = c('word', 'pvalue')
      # print(pathways)
      pdf(paste("../figures/accessory-essential-",prev_row[1],prev_row[2],prev_row[3],prev_row[4],prev_row[5],".pdf",sep=''))
      print(ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
              geom_abline(slope = 0, intercept = -log10(0.05), color="red")+
              # ggtitle("Core pathways") +
              theme(text = element_text(size=15)))
      dev.off()
    }
    start = ind
  }
  ind = ind + 1
}
genes = kegg[start:(ind-1),1]

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
    if (!is.na(item) & item == db[i,1])
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
  fisherresult = fisher.test(matrix(c(enrichment[i,1],enrichment[i,2],(len_genes-enrichment[i,1]), (len_db-enrichment[i,2])),nrow=2,byrow=TRUE), alternative = "greater")
  enrichment[i,3] = fisherresult$p.value
  
  # x = enrichment[i,1]
  # m = enrichment[i,2]
  # n = len_db - enrichment[i,2]
  # k = len_genes
  # enrichment[i,3] = phyper(x,m,n,k, lower.tail = TRUE)
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
if (nrow(pathways) > 0)
{
  colnames(pathways) = c('word', 'pvalue')
  # print(pathways)
  pdf(paste("../figures/accessory-essential-",prev_row[1],prev_row[2],prev_row[3],prev_row[4],prev_row[5],".pdf",sep=''))
  print(ggplot(data=pathways, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
          geom_abline(slope = 0, intercept = -log10(0.05), color="red")+
          # ggtitle("Core pathways") +
          theme(text = element_text(size=15)))
  dev.off()
}