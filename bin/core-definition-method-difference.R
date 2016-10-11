library(ggplot2)
codons <- "../results/word-enrichment/ancestral-intersection-info.txt"
codon_table <- read.table(codons, header = TRUE)
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "less")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment/ancestral-intersection-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-3,-9, -10,-11,-12,-14,-15,-18, -24,-25,-26,-29),]
names(df) = c("word", "pvalue")
pdf("../figures/ancestral-intersection-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Ancestral vs. intersection") +
  theme(text = element_text(size=15))
dev.off()


codons <- "../results/word-enrichment/dollo-ancestral-info.txt"
codon_table <- read.table(codons, header = TRUE)
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "less")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment/dollo-ancestral-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4,-7,-20,-21),]
names(df) = c("word", "pvalue")
pdf("../figures/dollo-ancestral-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Dollo vs. ancestral") +
  theme(text = element_text(size=15))
dev.off()


codons <- "../results/word-enrichment/dollo-intersection-info.txt"
codon_table <- read.table(codons, header = TRUE)
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "less")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment/dollo-intersection-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/dollo-intersection-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Dollo vs. intersection") +
  theme(text = element_text(size=15))
dev.off()