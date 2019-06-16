library("UpSetR")
library(stringr)
# listinput <- list(a=c(1,1,1,1,2,1,2,3,3,3), b=c(1,4,5), d=c(1,2,5,6,7), e=c(4,8,9))
# upset(fromList(listinput))

clusters_path <- "../results/hieranoid/clusters.txt"
# essentiality_path <- "../results/insertion-indices/normalised-insertion-indices/"
essentiality_path <- "../results/biases/dbscan/"
k12_path <- "../results/ecogene-k12.txt"
list_of_files <- list.files(path=essentiality_path, full.names=T, recursive=FALSE)
list_of_essential_genes = list()
for (filename in list_of_files)
{
  iis <- as.matrix(read.table(filename))
  list_of_essential_genes[strsplit(basename(filename),"\\.")[[1]][1]] <- list(iis[,1][iis[,3]=='essential'])
}
iis = as.matrix(read.table(k12_path))
list_of_essential_genes["b"] <- list(iis[,1])
presence = rev(list("BN373"=c(), "ERS227112"=c(), "ROD"=c(), "SL1344"=c(), "SL3261"=c(), "STMMW"=c(), "STM"=c(), "SEN"=c(),
                "t"=c(), "EC958"=c(), "NCTC13441"=c(), "BW25113"=c(), "b"=c()))
essentiality = rev(list("BN373"=c(), "ERS227112"=c(), "ROD"=c(), "SL1344"=c(), "SL3261"=c(), "STMMW"=c(), "STM"=c(), "SEN"=c(),
                "t"=c(), "EC958"=c(), "NCTC13441"=c(), "BW25113"=c(), "b"=c()))
cluster <- readLines(clusters_path)
numclusters = length(cluster)
for (clusterindex in (1:numclusters))
{
  line = cluster[clusterindex]
  genes = unlist(strsplit(line , "\t"))
  match_result = str_match_all(line, '([[:alnum:]]+_|[[:alpha:]]+)[[:alnum:]]+')
  species = gsub('_','',match_result[[1]][,2])
  # match_result = str_match_all(line, '[[,(]](([[:alnum:]]+)_[[:alnum:]]+):')
  # temp_result = str_match_all(line, '[[,(]](([[:alpha:]]+)[[:digit:]]+):')
  # match_result[[1]] <- rbind(match_result[[1]], temp_result[[1]])
  # match_result <- match_result[[1]]
  for (i in (1:length(species)))
  {
    name = species[i]
    if (name %in% names(presence)) 
    {
      presence[[name]] <- rbind(presence[[name]], clusterindex)
      gene = genes[i]
      if (gene %in% list_of_essential_genes[[name]])
      {
        essentiality[[name]] <- rbind(essentiality[[name]], clusterindex)
      }
    }
  }
}

names(presence) <- rev(c("Klebsiella pneumoniae Ecl8", "Klebsiella pneumoniae RH201207",
                     "Citrobacter rodentium ICC168", "Salmonella typhimurium SL1344", "Salmonella typhimurium SL3261",
                     "Salmonella typhimurium D23580", "Salmonella typhimurium A130", "Salmonella enteritidis P125109",
                     "Salmonella typhi Ty2", "Escherichia coli ST131 EC958", "Escherichia coli UPEC ST131"
                     , "Escherichia coli BW25113", "Escherichia coli K-12 MG1655"))
names(essentiality) <- rev(c("Klebsiella pneumoniae Ecl8", "Klebsiella pneumoniae RH201207",
                     "Citrobacter rodentium ICC168", "Salmonella typhimurium SL1344", "Salmonella typhimurium SL3261",
                     "Salmonella typhimurium D23580", "Salmonella typhimurium A130", "Salmonella enteritidis P125109",
                     "Salmonella typhi Ty2", "Escherichia coli ST131 EC958", "Escherichia coli UPEC ST131",
                     "Escherichia coli BW25113", "Escherichia coli K-12 MG1655"))

pdf("../figures/upsetr.pdf")
upset(fromList(presence), nsets = 14, order.by="freq", nintersects=25, keep.order = T, sets=names(presence))
      #, show.numbers = "no")
upset(fromList(essentiality), nsets = 14, order.by="freq", nintersects=25, keep.order = T, sets=names(essentiality))
dev.off()

# cs17ess=c()
# for (item in essentiality$`Escherichia coli ETEC CS17`)
# {
#   if (!(item %in% essentiality$`Klebsiella pneumoniae Ecl8` | item %in% essentiality$`Klebsiella pneumoniae RH201207` |
#         item %in% essentiality$`Enterobacter cloacae NCTC 9394` | item %in% essentiality$`Citrobacter rodentium ICC168` |
#         item %in% essentiality$`Salmonella Typhimurium SL1344` | item %in% essentiality$`Salmonella Typhimurium SL3261` |
#         item %in% essentiality$`Salmonella Typhimurium A130` | item %in% essentiality$`Salmonella Typhimurium D23580` |
#         item %in% essentiality$`Salmonella Enteritidis P125109` | item %in% essentiality$`Salmonella Typhi Ty2` |
#         item %in% essentiality$`Escherichia coli UPEC ST131` | item %in% essentiality$`Escherichia coli ETEC H10407` |
#         item %in% essentiality$`Escherichia coli K-12 MG1655`))
#   {
#     match_result = str_match(cluster[item], '(CS17_[[:alnum:]]+)')[[1]]
#     cs17ess = c(cs17ess,match_result)
#   }
# }
# library(Biostrings)
# cs17fasta = readDNAStringSet('../data/fasta-dna//chromosome/CS17.fasta')
# cs17towrite = c()
# for (item in cs17ess)
# {
#   cs17towrite = c(cs17towrite,cs17fasta@ranges@NAMES[startsWith(cs17fasta@ranges@NAMES,item)])
# }
# write.table(cs17towrite, file='../results/cs17ess.txt', quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)