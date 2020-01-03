library(stringr)
clusters_path <- "../results/merge-clust-plot"
list_of_files <- list.files(path=clusters_path, full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "SL3261", "STMMW", "t", "b", "BW25113", "EC958")
genuses = c(1, 2, 3, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 2, 2, 2)
names(genuses) <- names
numspecies = length(names)
file_II = list()
file_size = list()
file_group = list()
tbl_clst_species = list()

for (filename in list_of_files)
{
  clustspecies = c()
  cluster <- as.matrix(read.table(filename))
  i_sum = 0
  l_sum = 0
  cluster_size = nrow(cluster)
  clust_with_ii_size = 0
  for (i in (1:cluster_size))
  {
    match = str_match(cluster[i,2], "([[:graph:]]+)\\_[[:alnum:]]+")[2]
    if (is.na(match))
    {
      match = str_match(cluster[i,2], "([[:alpha:]]+)[[:digit:]]+")[2]
    }
    if (match %in% names)
    {
      clustspecies = c(clustspecies, match) 
    }
    i_sum = i_sum + as.numeric(cluster[i, 5])
    l_sum = l_sum + (as.numeric(cluster[i, 4]) - as.numeric(cluster[i, 3]) + 1) * 3
    clust_with_ii_size = clust_with_ii_size  + 1
  }
  if (l_sum > clust_with_ii_size * 60)
  {
    file_II[basename(filename)] = i_sum / clust_with_ii_size
    file_size[basename(filename)] = cluster_size
    unique_clustspecies = unique(clustspecies)
    one_or_more = length(unique_clustspecies)
    greater_than_one = length(table(clustspecies)[table(clustspecies)>1])
    cluster_genuses = c()
    for (item in unique_clustspecies)
    {
      cluster_genuses <- c(cluster_genuses, genuses[unique_clustspecies])
    }
    # if (as.numeric(one_or_more) <= 0.3 * numspecies)
    if (length(unique(cluster_genuses)) <= 1)
    {
      file_group[basename(filename)] = 'ORFan'
    }
    else if (as.numeric(one_or_more) - as.numeric(greater_than_one) >= 0.7 * as.numeric(one_or_more))
    {
      file_group[basename(filename)] = 'Single-copy'
    }
    else
    {
      file_group[basename(filename)] = 'Multiple-copy'
      if (min(table(clustspecies))==1)
      {
        tbl_clst_species[basename(filename)] <- list(table(clustspecies))
      }
    }
  }
}
insertion_index <- sapply(file_II, function(x){as.numeric(x[1])})
size_index <- sapply(file_size, function(x){as.numeric(x[1])})
group_index <- file_group

insertion_indices <- '/home/fatemeh/EnTrI/results/biases/dbscan'
list_of_files <- list.files(path=insertion_indices, full.names=T, recursive=FALSE)
essentials <- c()
for (filename in list_of_files)
{
  iis <- read.table(filename)
  temp <- as.character(iis[,1][iis[,3]=='essential'])
  essentials <- c(essentials, temp)
}

for (name in names(tbl_clst_species))
{
  filename = paste(clusters_path,name,sep='/')
  cluster <- as.character(read.table(filename)$V2)
  essentiality <- cluster %in% essentials
  if(sum(essentiality))
  {
    clustspecies <- c()
    for (item in cluster)
    {
      match = str_match(item, "([[:graph:]]+)\\_[[:alnum:]]+")[2]
      if (is.na(match))
      {
        match = str_match(item, "([[:alpha:]]+)[[:digit:]]+")[2]
      }
      if (match %in% names)
      {
        clustspecies = c(clustspecies, match) 
      }
    }
    s <- names(table(clustspecies))[table(clustspecies)==1]
    singlecopy <- clustspecies %in% s
    if (sum(singlecopy&essentiality)==sum(singlecopy))
    {
      print(name)
      # print(cluster)
      # print(essentiality)
      # print(singlecopy)
    }
  }
}