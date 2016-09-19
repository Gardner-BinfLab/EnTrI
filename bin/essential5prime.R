codons <- "../results/insertion-indices/essential5prime-info.out"
codon_table <- read.table(codons, header = TRUE)
fishertable =data.frame()
for (i in seq(1,(nrow(codon_table)-1)))
{
  fisherresult = fisher.test(matrix(c(codon_table[i,2],codon_table[i,3],codon_table[nrow(codon_table),2],codon_table[nrow(codon_table),3]),nrow=2,byrow=TRUE), alternative = "greater")
  fishertable[i,1] = codon_table[i,1]
  fishertable[i,2] = fisherresult$p.value
}
fishertable[,3] = p.adjust(fishertable[,2], method = "BH")
fishertable = fishertable[order(fishertable[,3], fishertable[,2]),]
write.table(fishertable, file = "../results/insertion-indices/essential5prime-pval.txt", quote = FALSE, col.names = c("codon", "p-value", "corrected-p-value"), row.names = FALSE, sep="\t")
