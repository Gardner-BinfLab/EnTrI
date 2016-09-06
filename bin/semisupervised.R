library(stringr)
library("upclass")
datapath <- "../results/insertion-indices/check-biases/"
essentialitypath <- "../results/insertion-indices/normalised-insertion-indices/"
k12_path <- "../results/ecogene-k12.txt"
clusters_path <- "../results/hieranoid/hieranoid-result.txt"
list_of_files <- list.files(path=datapath, full.names=T, recursive=FALSE)
datatable = data.frame()
counter = 0
for (filename in list_of_files)
{
  thistable <- read.table(filename, header = FALSE)
  thistable[,7] = counter
  counter = counter + 1
  datatable = rbind(datatable, thistable)
}
datatable <- datatable[-3]
row.names(datatable) <- datatable[,1]
datatable <- datatable[,-1]
names(datatable) <- c("ii", "distance", "gc", "length", "strain#")
for (coli in (1:ncol(datatable)))
{
  
  if (max(datatable[,coli]) > 1)
  {
    datatable[,coli]=datatable[,coli]/max(datatable[,coli])
  }
}
classes <- sapply(row.names(datatable),function(x) NULL)
iis=data.frame()
list_of_essential_genes = c()
species = c('b')
for (filename in list_of_files)
{
  iis <- as.matrix(read.table(filename))
  list_of_essential_genes <- c(list_of_essential_genes, as.vector(iis[,1][iis[,3]=='essential']))
  species = c(species, strsplit(basename(filename),"\\.")[[1]][1])
}
iis = as.matrix(read.table(k12_path))[,1]
list_of_essential_genes <- c(list_of_essential_genes,iis)

cluster <- readLines(clusters_path)
numclusters = length(cluster)
for (clusterindex in (1:numclusters))
{
  line = cluster[clusterindex]
  match_result = str_match_all(line, '[[,(]](([[:alnum:]]+)_[[:alnum:]]+):')
  temp_result = str_match_all(line, '[[,(]](([[:alpha:]]+)[[:digit:]]+):')
  match_result[[1]] <- rbind(match_result[[1]], temp_result[[1]])
  match_result <- match_result[[1]]
  if (setequal(match_result[,3], species) & length(match_result[,3]) == length(species))
  {
    sum = 0
    for (item in match_result[,2])
    {
      if (item %in% list_of_essential_genes)
      {
        sum = sum + 1
      }
    }
    if (sum == 0)
    {
      for (item in match_result[,2])
      {
        classes[item] = 0
      }
    }
    else if (sum == length(species))
    {
      for (item in match_result[,2])
      {
        classes[item] = 1
      }
    }
  }
}

xtrain = data.frame()
cltrain = c()
xtest = data.frame()
cltest = c()
for (i in (1:nrow(datatable)))
{
  if (!is.null(classes[[i]]))
  {
    xtrain = rbind(xtrain , datatable[i,])
    cltrain = c(cltrain, classes[[i]])
  }
  else
  {
    xtest = rbind(xtest , datatable[i,])
    cltest = c(cltest, classes[[i]])
  }
}
xtrain = as.matrix(xtrain)
xtest=as.matrix(xtest)
fitupmodels <- upclassify(xtrain, cltrain, xtest, cltest)
plot(fitupmodels)

for (i in (1:nrow(datatable)))
{
  if (is.null(classes[[i]]))
  {
    classes[i] <- fitupmodels$Best$test$cl[i]-1
  }
}

for (i in (1:nrow(datatable)))
{
  match_result = str_match(rownames(datatable)[i], '(([[:alnum:]]+)_[[:alnum:]]+):')
  if (is.na(match_result))
  {
    temp_result = str_match(rownames(datatable)[i], '(([[:alpha:]]+)[[:digit:]]+):')
  }
}
