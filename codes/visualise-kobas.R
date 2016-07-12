library(stringr)
library(ggplot2)
lineset <- readLines("../KOBAS/beneficial-losses/t.kobas")
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
df = df[-seq(21,nrow(df),1),]
pdf("../results/pathway-enrichment-ty2-beneficiallosses.pdf")
ggplot(data=df, aes(x=reorder(pathway,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Pathway',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Beneficial losses") +
  theme(text = element_text(size=15))
dev.off()