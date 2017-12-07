library(stringr)
library(dcGOR)
Pfam <- dcRDataLoader('Pfam')
background_ids = c()
domtblout_file = readLines("../results/dcgo/seqdb_Pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    background_ids = c(background_ids, match)
  }
}

pfam_ids = c()
domtblout_file = readLines("../results/dcgo/beneficial-loss_pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
benloss <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP", background=background_ids)
# pdf("../result/core-essential.pdf")
# visEnrichment(benloss, num_top_nodes = 20)
# dev.off()
a <- view(benloss, top_num=NA, sortBy="adjp", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$adjp)
colnames(b) = c("Go", "name", "adjusted_P-value")
c= b[b$`adjusted_P-value`<0.01,]
write.table(c, file="../results/dcgo/beneficial-loss.txt", row.names=FALSE, quote = FALSE)

pfam_ids = c()
domtblout_file = readLines("../results/dcgo/non-essential_pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
noness <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP", background=background_ids)
# pdf("../result/core-essential.pdf")
# visEnrichment(noness, num_top_nodes = 20)
# dev.off()
a <- view(noness, top_num=NA, sortBy="adjp", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$adjp)
colnames(b) = c("Go", "name", "adjusted_P-value")
c= b[b$`adjusted_P-value`<0.01,]
write.table(c, file="../results/dcgo/non-essential.txt", row.names=FALSE, quote = FALSE)

pfam_ids = c()
domtblout_file = readLines("../results/dcgo/essential_pfam-A.domtblout")
for (line in domtblout_file)
{
  match = str_match(line, "PF[:digit:]+")
  if(!is.na(match))
  {
    pfam_ids = c(pfam_ids, match)
  }
}
ess <- dcEnrichment(pfam_ids, domain="Pfam", ontology="GOBP", background=background_ids)
# pdf("../result/core-essential.pdf")
# visEnrichment(noness, num_top_nodes = 20)
# dev.off()
a <- view(ess, top_num=NA, sortBy="adjp", details=TRUE)
b = cbind.data.frame(a$term_id, a$term_name, a$adjp)
colnames(b) = c("Go", "name", "adjusted_P-value")
c= b[b$`adjusted_P-value`<0.01,]
write.table(c, file="../results/dcgo/essential.txt", row.names=FALSE, quote = FALSE)
