library("UpSetR")
library(stringr)
# listinput <- list(a=c(1,1,1,1,2,1,2,3,3,3), b=c(1,4,5), d=c(1,2,5,6,7), e=c(4,8,9))
# upset(fromList(listinput))

clusters_path <- "../results/merge-clust-plot"
list_of_files <- list.files(path=clusters_path, full.names=T, recursive=FALSE)
presence = list("ROD"=c(), "CS17"=c(), "ENC"=c(), "ETEC"=c(), "NCTC13441"=c(), "ERS227112"=c(), "BN373"=c(), "SEN"=c(), "STM"=c(),
                "SL1344"=c(), "STMMW"=c(), "t"=c(), "b"=c())
essentiality = list("ROD"=c(), "CS17"=c(), "ENC"=c(), "ETEC"=c(), "NCTC13441"=c(), "ERS227112"=c(), "BN373"=c(), "SEN"=c(), "STM"=c(),
                    "SL1344"=c(), "STMMW"=c(), "t"=c(), "b"=c())
list_of_essential_genes = list("ROD"=c(), "CS17"=c(), "ENC"=c(), "ETEC"=c(), "NCTC13441"=c(), "ERS227112"=c(), "BN373"=c(), "SEN"=c(), "STM"=c(),
                    "SL1344"=c(), "STMMW"=c(), "t"=c(), "b"=c())
for (filename in list_of_files)
{
  cluster <- as.matrix(read.table(filename))
  cluster_size = nrow(cluster)
  for (i in (1:cluster_size))
  {
    if (as.numeric(cluster[i,5]) >= 0)
    {
      match = str_match(cluster[i,2], "([[:graph:]]+)\\_[[:alnum:]]+")[2]
      if (is.na(match))
      {
        match = str_match(cluster[i,2], "([[:alpha:]]+)[[:digit:]]+")[2]
      }
      if (match %in% names(presence))
      {
        presence[[match]] = c(presence[[match]], strsplit(basename(filename),"\\.")[[1]][1])
        if (as.numeric(cluster[i,5]) < 0.2)
        {
          essentiality[[match]] = c(essentiality[[match]], strsplit(basename(filename),"\\.")[[1]][1])
          list_of_essential_genes[[match]] = c(list_of_essential_genes[[match]], cluster[i,2])
        }
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
upset(fromList(presence), nsets = 13, order.by="freq", nintersects=66)
upset(fromList(essentiality), nsets = 13, order.by="freq", nintersects=51)
