library(ComplexHeatmap)
toplot <- read.table('~/EnTrI/results/venn-entero-deg-endo/heatmap.tsv', header = TRUE)
rownames(toplot) = make.names(toplot$gene, unique=TRUE)
toplot <- toplot[,seq(2,6)]
toplot <- toplot[ order(toplot$Enterobacteriaceae, toplot$Endosymbiont, toplot$Gammaproteobacteria, toplot$Proteobacteria, toplot$Bacteria,
                        row.names(toplot), decreasing = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), method = 'radix'), ]
pdf("~/EnTrI/figures/essential_genes_heatmap.pdf", height = 60, width = 3)
Heatmap(toplot, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = 'black'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()

# dict <- c('Bacillus subtilis 168', 'Haemophilus influenzae Rd KW20', 'Helicobacter pylori 26695', 'Acinetobacter baylyi ADP1',
#           'Salmonella enterica serovar Typhi', 'Caulobacter crescentus', 'Burkholderia thailandensis E264', 'Sphingomonas wittichii RW1',
#           'Salmonella enterica serovar Typhimurium SL1344', 'Burkholderia pseudomallei', 'Streptococcus pyogenes NZ131',
#           'Rhodopseudomonas palustris CGA009', 'Agrobacterium fabrum str. C58', 'Escherichia coli ST131 strain EC958',
#           'Staphylococcus aureus N315', 'Mycoplasma genitalium G37', 'Salmonella typhimurium LT2', 'Mycoplasma pulmonis UAB CTIP',
#           'Staphylococcus aureus NCTC 8325', 'Streptococcus sanguinis', 'Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S',
#           'Shewanella oneidensis MR-1', 'Salmonella enterica serovar Typhi Ty2', 'Pseudomonas aeruginosa PAO1',
#           'Porphyromonas gingivalis ATCC 33277', 'Streptococcus agalactiae A909', 'Brevundimonas subvibrioides ATCC 15264',
#           'Vibrio cholerae N16961', 'Streptococcus pneumoniae', 'Francisella novicida U112', 'Pseudomonas aeruginosa UCBPP-PA14',
#           'Escherichia coli MG1655', 'Bacteroides thetaiotaomicron VPI-5482', 'Mycobacterium tuberculosis H37Rv',
#           'Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819', 'Bacteroides fragilis 638R', 'Streptococcus pyogenes MGAS5448',
#           'Synechococcus elongatus PCC 7942', 'Acinetobacter baumannii ATCC 17978', 'Bacillus thuringiensis BMB171')
# names(dict) <- c('DEG1001', 'DEG1005', 'DEG1008', 'DEG1013', 'DEG1016', 'DEG1020', 'DEG1024', 'DEG1028', 'DEG1032', 'DEG1035', 'DEG1038',
#                  'DEG1041', 'DEG1045', 'DEG1048', 'DEG1002', 'DEG1006', 'DEG1011', 'DEG1014', 'DEG1017', 'DEG1021', 'DEG1026', 'DEG1029',
#                  'DEG1033', 'DEG1036', 'DEG1039', 'DEG1042', 'DEG1046', 'DEG1003', 'DEG1007', 'DEG1012', 'DEG1015', 'DEG1019', 'DEG1023',
#                  'DEG1027', 'DEG1031', 'DEG1034', 'DEG1037', 'DEG1040', 'DEG1043', 'DEG1047')
dict <- c('Mycobacterium tuberculosis H37Rv', 'Synechococcus elongatus PCC 7942', 'Porphyromonas gingivalis ATCC 33277',
          'Bacteroides thetaiotaomicron VPI-5482', 'Bacteroides fragilis 638R', 'Haemophilus influenzae Rd KW20',
          'Escherichia coli ST131 strain EC958', 'Escherichia coli MG1655', 'Salmonella enterica serovar Typhi Ty2',
          'Salmonella enterica serovar Typhi', 'Salmonella enterica serovar Typhimurium SL1344', 'Salmonella typhimurium LT2',
          'Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S', 'Vibrio cholerae N16961', 'Shewanella oneidensis MR-1',
          'Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa UCBPP-PA14', 'Acinetobacter baumannii ATCC 17978', 'Acinetobacter baylyi ADP1',
          'Francisella novicida U112', 'Burkholderia thailandensis E264', 'Burkholderia pseudomallei K96243', 'Sphingomonas wittichii RW1',
          'Rhodopseudomonas palustris CGA009', 'Agrobacterium fabrum str. C58', 'Caulobacter crescentus',
          'Brevundimonas subvibrioides ATCC 15264', 'Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819', 'Helicobacter pylori 26695',
          'Mycoplasma genitalium G37', 'Mycoplasma pulmonis UAB CTIP', 'Staphylococcus aureus NCTC 8325', 'Staphylococcus aureus N315',
          'Bacillus subtilis 168', 'Bacillus thuringiensis BMB171', 'Streptococcus agalactiae A909', 'Streptococcus pyogenes NZ131',
          'Streptococcus pyogenes MGAS5448', 'Streptococcus sanguinis', 'Streptococcus pneumoniae')
names(dict) <- c('DEG1027', 'DEG1040', 'DEG1039', 'DEG1023', 'DEG1034', 'DEG1005', 'DEG1048', 'DEG1019', 'DEG1033', 'DEG1016', 'DEG1032',
                 'DEG1011', 'DEG1026', 'DEG1003', 'DEG1029', 'DEG1036', 'DEG1015', 'DEG1043', 'DEG1013', 'DEG1012', 'DEG1024', 'DEG1035',
                 'DEG1028', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1031', 'DEG1008', 'DEG1006', 'DEG1014', 'DEG1017', 'DEG1002',
                 'DEG1001', 'DEG1047', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1021', 'DEG1007')
eggnog <- read.table('~/EnTrI/results/deg/fasta-proteins/seqdb.fasta.emapper.annotations', sep = '\t', quote = '')
eggnogplot <- data.frame(matrix(0, nrow = nrow(toplot), ncol = length(dict)))
rownames(eggnogplot) <- row.names(toplot)
colnames(eggnogplot) <- dict
for (genename in row.names(eggnogplot))
{
  genomes = eggnog[,1][eggnog[,5]==toupper(genename)]
  for (item in genomes)
  {
    gene = substr(toString(item),1,7)
    eggnogplot[genename, dict[gene]] = 1
  }
}
pdf("~/EnTrI/figures/essential_genes_status_heatmap.pdf", height = 60, width = 25)
Heatmap(eggnogplot, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = 'black'), col=c('white', 'black'), name = 'Essentiality',
        heatmap_legend_param=c('Non-essential', 'Essential'))
dev.off()