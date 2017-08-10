eggnogspath <- '../results/luca-vs-ancestrallyessential/b.fasta.emapper.annotations'
tbl <- read.table(eggnogspath, sep='\t', quote = "")
eggnogs <- tbl[,c(1,5,12)]
colnames(eggnogs) <- c('gene_id', 'gene_name', 'gene_function')

keggpath <- '../results/KEGG_modules/escherichia_coli_K-12_MG1655.dat'
keggs <- read.table(keggpath, header = TRUE)
keggs <- keggs[,-2]

annots <- merge(eggnogs, keggs, by = 'gene_id', all.x = TRUE)
annots <- aggregate(annots[4], annots[-4], FUN = function(X) paste(unique(X), collapse="/ "))
sorted_annots <- annots[order(annots$gene_id),]

write.table(sorted_annots, '../results/luca-vs-ancestrallyessential/ancestrally-essential.tsv', quote = FALSE, sep = "\t", row.names = FALSE)

# combine the results with DEG output: ~/DEG.zip and save them in ancestrally-essential.tsv.backup

deg <- read.table('../results/luca-vs-ancestrallyessential/ancestrally-essential.tsv.backup', header = TRUE, quote = "", sep="\t")
unique_ess <- deg[deg$No_essential_homologs<6,]
sorted_unique_ess <- unique_ess[order(unique_ess$path_name),]
write.table(sorted_unique_ess, '../results/luca-vs-ancestrallyessential/ancestrally-essential-unique.tsv', quote = FALSE, sep = "\t", row.names = FALSE)
