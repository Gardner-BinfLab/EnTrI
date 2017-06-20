library(ggplot2)
library(stringr)
codons <- "../results/word-enrichment2/beneficialloss-info.txt"
# db <- "../results/eggnog-mapper/gproNOG.annotations.tsv"
codon_table <- read.table(codons, header = TRUE, sep = '\t', quote="")
# eggnogdf <- read.table(db, header = FALSE, sep='\t', quote = "")
# eggnogs <- eggnogdf[,6]
# names(eggnogs) <- eggnogdf[,2]
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "greater")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment2/beneficialloss-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

# df = fishertable[-seq(25,nrow(fishertable),1),]
df = fishertable#[!is.na(fishertable$V1),]
df[,3] = -log10(df[,3])
df$V1 <- str_match(df$V1, '([[:graph:]]+([[:blank:]]+[[:graph:]]+){0,6})')[,1]
df = df[seq(1,30),-2]
# df = df[c(-2,-5,-9,-13),]
names(df) = c("word", "pvalue")
df$word[1] <- 'Not classified'
df$word[2] <- 'pct identical to transposase Yersinia pestis'
df$word[6] <- 'Inhibits RpoS proteolysis'
df$word[20] <- 'Influences nickel and cobalt homeostasis'
# df$word[30] <- 'Catalyzes the dephosphorylation of heptose(II)'
pdf("../figures/beneficialloss-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Group',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Beneficial losses")+
  theme(text = element_text(size=15))
dev.off()

codons <- "../results/word-enrichment2/essential-info.txt"
codon_table <- read.table(codons, header = TRUE, sep = '\t', quote="")
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "greater")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment2/essential-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

# df = fishertable[-seq(30,nrow(fishertable),1),]
df = fishertable[!is.na(fishertable$V1),]
df[,3] = -log10(df[,3])
df$V1 <- str_match(df$V1, '([[:graph:]]+([[:blank:]]+[[:graph:]]+){0,6})')[,1]
df = df[seq(1,30),-2]
#df = df[c(-4,-7,-20,-21),]
names(df) = c("word", "pvalue")
pdf("../figures/essential-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ggtitle("Essential genes") #+
  #theme(text = element_text(size=24)) +theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))
dev.off()


codons <- "../results/word-enrichment2/non-essential-info.txt"
codon_table <- read.table(codons, header = TRUE, sep = '\t', quote="")
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],(codon_table[nrow(codon_table),2]-codon_table[i,2]),
                                      (codon_table[nrow(codon_table),3])-codon_table[i,3]),nrow=2,byrow=TRUE), alternative = "greater")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/word-enrichment2/non-essential-pval.txt", quote = FALSE, col.names = c("word", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")

# df = fishertable[-seq(30,nrow(fishertable),1),]
df = fishertable[!is.na(fishertable$V1),]
df[,3] = -log10(df[,3])
df$V1 <- str_match(df$V1, '([[:graph:]]+([[:blank:]]+[[:graph:]]+){0,6})')[,1]
df = df[seq(1,30),-2]
#df = df[c(-4),]
names(df) = c("word", "pvalue")
pdf("../figures/non-essential-pval.pdf")
ggplot(data=df, aes(x=reorder(word,pvalue),y=pvalue))+geom_bar(stat="identity")+ coord_flip()+labs(x='Word',y='-log10(P-value)') +
  geom_abline(slope = 0, intercept = -log10(0.05), color="red")+ylim(c(0,max(df[1,2], -log10(0.05))))+ggtitle("Non-essential genes") #+
  #theme(text = element_text(size=24)) +theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))
dev.off()