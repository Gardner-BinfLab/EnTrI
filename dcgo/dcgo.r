library(stringr)
library(dcGOR)
pfam_ids = c()
domtblout_file = readLines("core-essential_Pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
Pfam <- dcRDataLoader('Pfam')
core_essential_out <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP")
#pdf("res/core-essential.pdf")
visEnrichment(core_essential_out, num_top_nodes = 20)
#dev.off()
a <- view(core_essential_out, top_num=20, sortBy="pvalue", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$pvalue)
colnames(b) = c("Go", "name", "P-value")
write.table(b, file="res/core-essential.txt", row.names=FALSE, quote = FALSE)

pfam_ids = c()
domtblout_file = readLines("core-never-essential_Pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
core_essential_out <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP")
#pdf("res/core-never-essential.pdf")
visEnrichment(core_essential_out, num_top_nodes = 20)
#dev.off()
a <- view(core_essential_out, top_num=20, sortBy="pvalue", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$pvalue)
colnames(b) = c("Go", "name", "P-value")
write.table(b, file="res/core-never-essential.txt", row.names=FALSE, quote = FALSE)

pfam_ids = c()
domtblout_file = readLines("accessory-essential_Pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
core_essential_out <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP")
#pdf("res/accessory-essential.pdf")
visEnrichment(core_essential_out, num_top_nodes = 20)
#dev.off()
a <- view(core_essential_out, top_num=20, sortBy="pvalue", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$pvalue)
colnames(b) = c("Go", "name", "P-value")
write.table(b, file="res/accessory-essential.txt", row.names=FALSE, quote = FALSE)

pfam_ids = c()
domtblout_file = readLines("accessory-never-essential_Pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
core_essential_out <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP")
#pdf("res/accessory-never-essential.pdf")
visEnrichment(core_essential_out, num_top_nodes = 20)
#dev.off()
a <- view(core_essential_out, top_num=20, sortBy="pvalue", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$pvalue)
colnames(b) = c("Go", "name", "P-value")
write.table(b, file="res/accessory-never-essential.txt", row.names=FALSE, quote = FALSE)