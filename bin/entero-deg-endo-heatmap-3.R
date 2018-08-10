library(ComplexHeatmap)
library(stringr)
library(UpSetR)
toplot <- read.table('~/EnTrI/results/venn-entero-deg-endo/heatmap.tsv', header = TRUE)
rownames(toplot) = make.names(toplot$gene, unique=TRUE)
toplot <- toplot[,seq(2,6)]
toplot <- toplot[ order(toplot$Enterobacteriaceae, toplot$Endosymbiont, toplot$Gammaproteobacteria, toplot$Proteobacteria, toplot$Bacteria,
                        row.names(toplot), decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), method = 'radix'), ]
pdf("~/EnTrI/figures/essential_genes_heatmap.pdf", height = 60, width = 3)
Heatmap(toplot, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = 'gray'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()

upsetrlist <- list("Enterobacteriaceae"=c(), "Endosymbiont"=c(), "Gammaproteobacteria"=c(), "Proteobacteria"=c(), "Bacteria"=c())
for (gene in row.names(toplot))
{
  for (cl in names(upsetrlist))
  if(toplot[gene,cl])
    upsetrlist[[cl]]=c(upsetrlist[[cl]], gene)
}

pdf("~/EnTrI/figures/essential_genes_upsetr.pdf")
upset(fromList(upsetrlist), nsets = 5, order.by="freq", keep.order = T, sets=names(upsetrlist))
dev.off()

dict <- c('Mycobacterium tuberculosis H37Rv', 'Synechococcus elongatus PCC 7942', 'Porphyromonas gingivalis ATCC 33277',
          'Bacteroides thetaiotaomicron VPI-5482','Bacteroides fragilis 638R', 'Caulobacter crescentus',
          'Brevundimonas subvibrioides ATCC 15264', 'Rhodopseudomonas palustris CGA009', 'Agrobacterium fabrum str. C58',
          'Burkholderia thailandensis E264', 'Burkholderia pseudomallei K96243', 'Francisella novicida U112',
          'Shewanella oneidensis MR-1', 'Vibrio cholerae N16961', 'Haemophilus influenzae Rd KW20',
          'Escherichia coli ST131 strain EC958', 'Escherichia coli MG1655', 'Salmonella enterica serovar Typhimurium SL1344',
          'Salmonella enterica serovar Typhi Ty2', 'Salmonella enterica serovar Typhi', 'Acinetobacter baumannii ATCC 17978',
          'Pseudomonas aeruginosa PAO1', 'Helicobacter pylori 26695', 'Mycoplasma genitalium G37',
          'Mycoplasma pulmonis UAB CTIP', 'Streptococcus agalactiae A909', 'Streptococcus pyogenes MGAS5448',
          'Streptococcus pyogenes NZ131', 'Staphylococcus aureus NCTC 8325', 'Staphylococcus aureus N315',
          'Bacillus subtilis 168')
names(dict) <- c('DEG1027', 'DEG1040', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1020', 'DEG1046', 'DEG1041',
                 'DEG1045', 'DEG1024', 'DEG1035', 'DEG1012', 'DEG1029', 'DEG1003', 'DEG1005', 'DEG1048',
                 'DEG1019', 'DEG1032', 'DEG1033', 'DEG1016', 'DEG1043', 'DEG1036', 'DEG1008', 'DEG1006',
                 'DEG1014', 'DEG1042', 'DEG1037', 'DEG1038', 'DEG1017', 'DEG1002', 'DEG1001')
eggnog <- read.table('~/EnTrI/results/deg/fasta-proteins/seqdb.fasta.emapper.annotations', sep = '\t', quote = '')
eggnogplot <- data.frame(matrix(0, nrow = nrow(toplot), ncol = length(dict)))
rownames(eggnogplot) <- row.names(toplot)
colnames(eggnogplot) <- dict
for (genename in row.names(eggnogplot))
{
  genomes = eggnog[,1][eggnog[,5]==toupper(genename)]
  for (item in genomes)
  {
    if (startsWith(item, 'exDEG'))
    {
      gene = substr(toString(item),3,9)
      eggnogplot[genename, dict[gene]] = 2
    }
    else
    {
      gene = substr(toString(item),1,7)
      eggnogplot[genename, dict[gene]] = 1
    }
  }
}
write.table(eggnogplot, file="~/EnTrI/results/venn-entero-deg-endo/deg-heatmap.tsv", quote = FALSE, sep='\t')
pdf("~/EnTrI/figures/essential_genes_status_deg_heatmap.pdf", height = 65, width = 25)
Heatmap(eggnogplot, cluster_rows = FALSE, cluster_columns = FALSE, column_names_max_height=unit(13, "cm"), rect_gp = gpar(col = 'gray'), col=c('white', 'black', 'gray'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential', 'Not analysed'))
dev.off()

dict <- c('Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
          'Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans', 'Sodalis glossinidius str. morsitans',
          'Sodalis praecaptivus strain HS1', 'Candidatus Sodalis pierantonius str. SOPE', 'Secondary endosymbiont of Ctenarytaina eucalypti',
          'Candidatus Moranella endobia PCIT', 'Secondary endosymbiont of Heteropsylla cubana',
          'Candidatus Baumannia cicadellinicola strain BGSS', 'Candidatus Baumannia cicadellinicola strain B-GSS', 
          'Baumannia cicadellinicola str. Hc (Homalodisca coagulata)', 'Candidatus Blochmannia chromaiodes str. 640',
          'Candidatus Blochmannia pennsylvanicus str. BPEN', 'Candidatus Blochmannia floridanus', 'Candidatus Blochmannia vafer str. BVAF',
          'Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675',
          'Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757', 'Buchnera aphidicola str. Bp (Baizongia pistaciae)',
          'Buchnera aphidicola BCc', 'Buchnera aphidicola (Cinara tujafilina)', 'Buchnera aphidicola str. F009 (Myzus persicae)',
          'Buchnera aphidicola str. W106 (Myzus persicae)', 'Buchnera aphidicola str. USDA (Myzus persicae)',
          'Buchnera aphidicola str. G002 (Myzus persicae)', 'Buchnera aphidicola str. Ua (Uroleucon ambrosiae)',
          'Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)', 'Buchnera aphidicola str. APS (Acyrthosiphon pisum)',
          'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)',
          'Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)',
          'Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)', 'Buchnera aphidicola str. 5A (Acyrthosiphon pisum)',
          'Buchnera aphidicola str. Sg (Schizaphis graminum)')
names(dict) <- c('BA000021', 'WIGMOR', 'SG', 'Sant', 'SOPEG', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BCHRO640', 'BPEN', 'Bfl',
                 'BVAF', 'BTURN675', 'BOBLI757', 'bbp', 'BCc', 'BCTU', 'BUMPF009', 'BUMPW106', 'BUMPUSDA', 'BUMPG002', 'BUAMB', 'BAKON',
                 'BA000003', 'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUsg')
eggnog <- read.table('~/EnTrI/results/endosymbionts/fasta-proteins/seqdb.fasta.emapper.annotations', sep = '\t', quote = '')
eggnogplot <- data.frame(matrix(0, nrow = nrow(toplot), ncol = length(dict)))
rownames(eggnogplot) <- row.names(toplot)
colnames(eggnogplot) <- dict
for (genename in row.names(eggnogplot))
{
  genomes = eggnog[,1][eggnog[,5]==toupper(genename)]
  for (item in genomes)
  {
    gene = str_match(item, "([[:graph:]]+?)\\_[[:graph:]]+")[2]
    if (is.na(gene))
    {
      gene = str_match(item, "([[:alpha:]]+)[[:digit:]]+")[2]
    }
    eggnogplot[genename, dict[gene]] = 1
  }
}
write.table(eggnogplot, file="~/EnTrI/results/venn-entero-deg-endo/endosymbionts-heatmap.tsv", quote = FALSE, sep='\t')
pdf("~/EnTrI/figures/essential_genes_status_endo_heatmap.pdf", height = 65, width = 25)
Heatmap(eggnogplot, cluster_rows = FALSE, cluster_columns = FALSE, column_names_max_height=unit(20, "cm"), rect_gp = gpar(col = 'gray'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()