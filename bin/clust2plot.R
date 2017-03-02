library(stringr)
# args <- commandArgs(trailingOnly = TRUE)
# clusters <- args[1]
clusters_path <- "../results/merge-clust-plot"
#clusters_path <- c("../results/merge-clust-plot", "../results/merge-clust-plot-without-ends/")
cutoff = 1.644854
for (cpitem in clusters_path)
{
  list_of_files <- list.files(path=cpitem, full.names=T, recursive=FALSE)
  names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "b", "BW25113", "EC958")
  genuses = c(1, 2, 3, 2, 2, 3, 3, 4, 4, 4, 4, 4, 2, 2, 2)
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
      }
    }
  }
  insertion_index <- sapply(file_II, function(x){as.numeric(x[1])})
  size_index <- sapply(file_size, function(x){as.numeric(x[1])})
  group_index <- file_group
  
  insertion_index=(insertion_index-mean(insertion_index))/sd(insertion_index-mean(insertion_index))
  
  # if (cpitem == clusters_path[1])
  #   pdf("../results/cluster-essentiality.pdf")
  # else
  #   pdf("../results/cluster-essentiality-without-ends.pdf")
  pdf("../figures/cluster-essentiality.pdf")
  
  m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
  temp <- split.screen(m)
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  
  h <- hist(insertion_index, breaks =seq(min(insertion_index),max(insertion_index)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,-cutoff, cutoff, Inf))
  screen(1)
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "NPEQ", main ="All gene classes", cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(-3,3), ylim=c(0,200), lty= "blank", axes=FALSE)
  axis(1, at=seq(-3,3,1), cex.axis=1.5)
  axis(2, at=seq(0,200,50), labels=c(0,NA,NA,NA,200), cex.axis=1.5)
  text(1,190, paste("n =", length(insertion_index)), lty=1, lwd=4, cex=1.15, bty="n")
  legend(1,180, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
         col=c("darkmagenta","turquoise4", "darkgoldenrod4"), bty="n")
  #lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  #lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
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
      if (insertion_index[item] < -cutoff)
        orfans_ben = c(orfans_ben, insertion_index[item])
      else if (insertion_index[item] < cutoff)
        orfans_non = c(orfans_non, insertion_index[item])
      else
        orfans_es = c(orfans_es, insertion_index[item])
    } 
    else if (group_index[item] == 'Single-copy')
    {
      single_occurrence = c(single_occurrence, insertion_index[item])
      if (insertion_index[item] < -cutoff)
        single_ben = c(single_ben, insertion_index[item])
      else if (insertion_index[item] < cutoff)
        single_non = c(single_non, insertion_index[item])
      else
        single_es = c(single_es, insertion_index[item])
    }
    else
    {
      multiple_copies = c(multiple_copies, insertion_index[item])
      if (insertion_index[item] < -cutoff)
        multiple_ben = c(multiple_ben, insertion_index[item])
      else if (insertion_index[item] < cutoff)
        multiple_non = c(multiple_non, insertion_index[item])
      else
        multiple_es = c(multiple_es, insertion_index[item])
    }
  }
  
  h <- hist(orfans, breaks =seq(min(insertion_index),(max(insertion_index)+1),0.02), plot = FALSE)
  cuts <- cut(h$breaks, c(-Inf,-cutoff, cutoff, Inf))
  screen(2)
  par(mar=c(5.1,2.5,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab=NA, ylab=NA, main ="Genus specific", cex.axis=1.5, cex.main = 1.5,
       xlim=c(-3,3), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(-2,2,2), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(0,NA,100), cex.axis=1.5)
  text(1,80, paste("n =", length(orfans)), lty=1, lwd=4, cex=1.15, bty="n")
  #lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  #lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  h <- hist(single_occurrence, breaks =seq(min(insertion_index),max(insertion_index)+1,0.02), plot = FALSE)
  cuts <- cut(h$breaks, c(-Inf,-cutoff, cutoff, Inf))
  screen(3)
  par(mar=c(5.1,1,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab=NA, ylab=NA, main ="Single copy", cex.axis=1.5, cex.main = 1.5,
       xlim=c(-3,3), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(-2,2,2), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(NA,NA,NA), cex.axis=1.5)
  text(1,80, paste("n =", length(single_occurrence)), lty=1, lwd=4, cex=1.15, bty="n")
  #lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  #lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  h <- hist(multiple_copies, breaks =seq(min(insertion_index),max(insertion_index)+1,0.02), plot = FALSE)
  cuts <- cut(h$breaks, c(-Inf,-cutoff, cutoff, Inf))
  screen(4)
  #par(mar=c(2,1,2,1))
  par(mar=c(5.1,1,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA, ylab=NA, main ="Multi-copy", cex.axis = 1.5, cex.main = 1.5,
       xlim=c(-3,3), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(-2,2,2), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(NA,NA,NA), cex.axis=1.5)
  text(1,80, paste("n =", length(multiple_copies)), lty=1, lwd=4, cex=1.15, bty="n")
  #lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  #lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  close.screen(all.screens = TRUE)
  
  dev.off()
}

print("Essential Orfans:")
names(file_group)[file_group=="ORFan" & file_II > cutoff & file_size > 1]
print("Essential Multi-copies:")
names(file_group)[file_group=="Multiple-copy" & file_II > cutoff]