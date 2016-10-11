library(ggplot2)
codons <- "../results/branch-root/ecoli-info.txt"
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
write.table(fishertable, file = "../results/branch-root/ecoli-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-3,-9, -10,-11,-12,-14,-15,-18, -24,-25,-26,-29),]
names(df) = c("word", "pvalue")
pdf("../figures/ecoli-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Ecoli") +
  theme(text = element_text(size=15))
dev.off()


codons <- "../results/branch-root/salmonella-info.txt"
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
write.table(fishertable, file = "../results/branch-root/salmonella-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4,-7,-20,-21),]
names(df) = c("word", "pvalue")
pdf("../figures/salmonella-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Salmonella") +
  theme(text = element_text(size=15))
dev.off()


codons <- "../results/branch-root/klebsiella-info.txt"
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
write.table(fishertable, file = "../results/branch-root/klebsiella-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/klebsiella-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Klebsiella") +
  theme(text = element_text(size=15))
dev.off()


codons <- "../results/branch-root/klebsiellaenterobacter-info.txt"
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
write.table(fishertable, file = "../results/branch-root/klebsiellaenterobacter-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/klebsiellaenterobacter-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("KlebsiellaEnterobacter") +
  theme(text = element_text(size=15))
dev.off()

codons <- "../results/branch-root/salmonellacitrobacter-info.txt"
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
write.table(fishertable, file = "../results/branch-root/salmonellacitrobacter-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/salmonellacitrobacter-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("SalmonellaCitrobacter") +
  theme(text = element_text(size=15))
dev.off()

codons <- "../results/branch-root/salmonellaecolicitrobacter-info.txt"
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
write.table(fishertable, file = "../results/branch-root/salmonellaecolicitrobacter-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

df = fishertable[-seq(33,nrow(fishertable),1),]
df[,3] = -log10(df[,3])
df = df[,-2]
# df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/salmonellaecolicitrobacter-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("SalmonellaEcoliCitrobacter") +
  theme(text = element_text(size=15))
dev.off()