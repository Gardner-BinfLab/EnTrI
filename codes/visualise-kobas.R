library(stringr)
library(ggplot2)
lineset <- readLines("../KOBAS/accessory-essential/t.kobas")
kegg_table <- c()
for (item in lineset)
  {
    if (!(startsWith(item, '#') | startsWith(item, '-') | item ==''))
    {
      kegg_table <- rbind(kegg_table, strsplit(item, '\t'))
    }
}
names <- c()
pval <- c()
for (item in kegg_table)
  {
    names <- c(names, item[1])
    pval <- c(pval,as.numeric(item[7]))
  }
df <- data.frame(pathway=names,pvalue=-log10(pval))
pdf("../results/pathway-enrichment-accessory-essential.pdf")
ggplot(data=df, aes(x=reorder(pathway,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='P-value') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")
dev.off()