library("UpSetR")
library(stringr)
# listinput <- list(a=c(1,1,1,1,2,1,2,3,3,3), b=c(1,4,5), d=c(1,2,5,6,7), e=c(4,8,9))
# upset(fromList(listinput))

clusters_path <- "../results/hieranoid-result.txt"
essentiality_path <- "../results/insertion-indices/normalised-insertion-indices/"
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
presence = list("ROD"=c(), "CS17"=c(), "ENC"=c(), "ETEC"=c(), "NCTC13441"=c(), "ERS227112"=c(), "BN373"=c(), "SEN"=c(), "STM"=c(),
                "SL1344"=c(), "STMMW"=c(), "t"=c(), "b"=c())
essentiality = list("ROD"=c(), "CS17"=c(), "ENC"=c(), "ETEC"=c(), "NCTC13441"=c(), "ERS227112"=c(), "BN373"=c(), "SEN"=c(), "STM"=c(),
                    "SL1344"=c(), "STMMW"=c(), "t"=c(), "b"=c())
cluster <- readLines(clusters_path)
numclusters = length(cluster)
for (clusterindex in (1:numclusters))
{
  line = cluster[clusterindex]
  match_result = str_match_all(line, '[[,(]](([[:alnum:]]+)_[[:alnum:]]+):')
  temp_result = str_match_all(line, '[[,(]](([[:alpha:]]+)[[:digit:]]+):')
  match_result[[1]] <- rbind(match_result[[1]], temp_result[[1]])
  match_result <- match_result[[1]]
  for (i in (1:nrow(match_result)))
  {
    name = match_result[i,3]
    if (name %in% names(presence)) 
    {
      presence[[name]] <- rbind(presence[[name]], clusterindex)
      gene = match_result[i,2]
      if (gene %in% list_of_essential_genes[[name]])
      {
        essentiality[[name]] <- rbind(essentiality[[name]], clusterindex)
      }
    }
  }
}

names(presence) <- c("Citrobacter", "EC_ETEC_CS17", "Enterobacter", "EC_ETEC_H10407", "EC_UPEC",
                     "KP_RH201207", "KP_Ecl8", "S_enteritidis", "S_typhimurium_A130",
                     "S_typhimurium_SL1344", "S_typhimurium_D23580", "S_typhi", "EC_k12")
names(essentiality) <- c("Citrobacter", "EC_ETEC_CS17", "Enterobacter", "EC_ETEC_H10407", "EC_UPEC",
                         "KP_RH201207", "KP_Ecl8", "S_enteritidis", "S_typhimurium_A130",
                         "S_typhimurium_SL1344", "S_typhimurium_D23580", "S_typhi", "EC_k12")
upset(fromList(presence), nsets = 13, order.by="freq", nintersects=20
      #, show.numbers = "no"
     )
upset(fromList(essentiality), nsets = 13, order.by="freq", nintersects=20)
