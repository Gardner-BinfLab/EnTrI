library(stringr)
library("MASS")
require(reshape2)  # this is the library that lets you flatten out data
require(ggplot2)
# args <- commandArgs(trailingOnly = TRUE)
# clusters <- args[1]
clusters_path <- "../results/merge-clust-plot"

list_of_files <- list.files(path=clusters_path, full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261")
genuses = c(1, 2, 3, 2, 2, 4, 4, 5, 5, 5, 5, 5, 5)
names(genuses) <- names
numspecies = length(names)
file_II = list()
file_size = list()
file_group = list()
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
    if (as.numeric(cluster[i,5]) >= 0)
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
    }
  }
}
insertion_index <- sapply(file_II, function(x){as.numeric(x[1])})
size_index <- sapply(file_size, function(x){as.numeric(x[1])})
group_index <- file_group

pdf("../figures/cluster-essentiality.pdf")

m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
# temp <- split.screen(m)

# mar.default <- c(5,4,4,2) + 0.1
# par(mar = mar.default + c(0, 1, 0, 0))

# screen(1)
ii = insertion_index
nG = length(ii)

#identify second maxima
h <- hist(ii, breaks=0:(max(ii)*50+1)/50,plot=FALSE)
maxindex <- which.max(h$density[3:length(h$density)])
maxval <- h$mids[maxindex+2]

#find inter-mode minima with loess
r <- floor(maxval *1000)
I = ii < r / 1000
h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density

m1 = h1$mids[which.min(predict(lo))]
m2 = h$mids[max(which(h$counts>5))]
I1 = ((ii < m1)&(ii > 0))
I2 = ((ii >= m1)&(ii < m2))

f1 = (sum(I1) + sum(ii == 0))/nG
f2 = (sum(I2))/nG

d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves
h <- hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", plot=FALSE)
lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
essen <- lower/10000
ambig <- upper/10000
noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])

orfans = c()
orfans_non = c()
orfans_es = c()
orfans_ben = c()
single_occurrence = c()
single_non = c()
single_es = c()
single_ben = c()
multiple_copies = c()
multiple_non = c()
multiple_es = c()
multiple_ben = c()
for (item in names(size_index))
{
  if (group_index[item] == 'ORFan')
  {
    orfans = c(orfans, insertion_index[item])
    if (insertion_index[item] < 0.2)
      orfans_es = c(orfans_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      orfans_non = c(orfans_non, insertion_index[item])
    else
      orfans_ben = c(orfans_ben, insertion_index[item])
  } 
  else if (group_index[item] == 'Single-copy')
  {
    single_occurrence = c(single_occurrence, insertion_index[item])
    if (insertion_index[item] < 0.2)
      single_es = c(single_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      single_non = c(single_non, insertion_index[item])
    else
      single_ben = c(single_ben, insertion_index[item])
  }
  else
  {
    multiple_copies = c(multiple_copies, insertion_index[item])
    if (insertion_index[item] < 0.2)
      multiple_es = c(multiple_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      multiple_non = c(multiple_non, insertion_index[item])
    else
      multiple_ben = c(multiple_ben, insertion_index[item])
  }
}

bucket<-list(orfans=orfans, single_occurrence=single_occurrence, multiple_copies=multiple_copies)
mlt <- melt(bucket)
for (i in seq(1,nrow(mlt)))
{
  if (mlt[i,2] == "orfans")
  {
    if (mlt[i,1] < essen)
    {
      mlt[i,2] = 1
    }
    else if(mlt[i,1] < ambig)
    {
      mlt[i,2] = 2
    }
    else if(mlt[i,1] < noness)
    {
      mlt[i,2] = 3
    }
    else
    {
      mlt[i,2] = 4
    }
  }
  else if(mlt[i,2] == "single_occurrence")
  {
    if (mlt[i,1] < essen)
    {
      mlt[i,2] = 5
    }
    else if(mlt[i,1] < ambig)
    {
      mlt[i,2] = 6
    }
    else if(mlt[i,1] < noness)
    {
      mlt[i,2] = 7
    }
    else
    {
      mlt[i,2] = 8
    }
  }
  else if(mlt[i,2] == "multiple_copies")
  {
    if (mlt[i,1] < essen)
    {
      mlt[i,2] = 9
    }
    else if(mlt[i,1] < ambig)
    {
      mlt[i,2] = 10
    }
    else if(mlt[i,1] < noness)
    {
      mlt[i,2] = 11
    }
    else
    {
      mlt[i,2] = 12
    }
  }
}
ggplot(mlt, aes(x=value)) + 
  geom_histogram(data=subset(mlt,L1==1 | L1==5 | L1==9), position = "stack", binwidth=0.05, fill='darkgoldenrod4', aes(alpha='a')) +
  geom_histogram(data=subset(mlt,L1==2 | L1==6 | L1==10), position = "stack", binwidth=0.05, fill='darkred', alpha=0.2) +
  geom_histogram(data=subset(mlt,L1==3 | L1==7 | L1==11), position = "stack", binwidth=0.05, fill='darkcyan', alpha=0.2) +
  geom_histogram(data=subset(mlt,L1==4 | L1==8 | L1==12), position = "stack", binwidth=0.05, fill='darkorchid4', alpha=0.2) +
  geom_histogram(data=subset(mlt,L1==5 | L1==9), position = "stack", binwidth=0.05, fill='darkgoldenrod4', aes(alpha='b')) +
  geom_histogram(data=subset(mlt,L1==6 | L1==10), position = "stack", binwidth=0.05, fill='darkred', alpha=0.6) +
  geom_histogram(data=subset(mlt,L1==7 | L1==11), position = "stack", binwidth=0.05, fill='darkcyan', alpha=0.6) +
  geom_histogram(data=subset(mlt,L1==8 | L1==12), position = "stack", binwidth=0.05, fill='darkorchid4', alpha=0.6) +
  geom_histogram(data=subset(mlt,L1==9), position = "stack", binwidth=0.05, aes(fill='00', alpha='c')) +
  geom_histogram(data=subset(mlt,L1==10), position = "stack", binwidth=0.05, aes(fill='01'), alpha=1) +
  geom_histogram(data=subset(mlt,L1==11), position = "stack", binwidth=0.05, aes(fill='02'), alpha=1) +
  geom_histogram(data=subset(mlt,L1==12), position = "stack", binwidth=0.05, aes(fill='03'), alpha=1) +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25,face="bold"), plot.title = element_text(face="bold", size=32)) +
  xlab("Insertion Index") + ylab("Frequency") +
  ggtitle("Clusters Insertion Index") +
  scale_fill_manual(name = 'Essentiality',  values =c('00'='darkgoldenrod4','01'='darkred','02'='darkcyan', '03'='darkorchid4'),
                    labels = c('Essential', 'Ambiguous', 'Non-essential', 'Beneficial loss')) +
  scale_alpha_manual(name = 'Conservation',  values =c('a'=0.2,'b'=0.6,'c'=1),
                    labels = c('Genus-specific', 'Single-copy', 'Multi-copy'))





dev.off()