library(stringr)
# args <- commandArgs(trailingOnly = TRUE)
# clusters <- args[1]
clusters_path <- c("../results/merge-clust-plot-without-ends/")
for (cpitem in clusters_path)
{
  list_of_files <- list.files(path=cpitem, full.names=T, recursive=FALSE)
  names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "b")
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
      one_or_more = length(unique(clustspecies))
      greater_than_one = length(table(clustspecies)[table(clustspecies)>1])
      if (as.numeric(one_or_more) <= 0.3 * numspecies)
      {
        file_group[basename(filename)] = 'ORFan'
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
  
  #pdf("../results/cluster-essentiality-without-ends-conserved-unconserved.pdf")
  
  m <- rbind(c(0,1,0.5,1), c(0, 0.5, 0, 0.5), c(0.5, 1, 0, 0.5))
  temp <- split.screen(m)
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  
  h <- hist(insertion_index, breaks =seq(min(insertion_index),max(insertion_index)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(1)
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "Insertion Index", main ="All clusters", cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,200,50), labels=c(0,NA,NA,NA,200), cex.axis=1.5)
  text(2.85,190, paste("n =", length(insertion_index)), lty=1, lwd=4, cex=1.15, bty="n")
  legend(2.4,180, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
         col=c("darkgoldenrod4","turquoise4", "darkmagenta"), bty="n")
  lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  orfans = c()
  orfans_non = c()
  orfans_es = c()
  orfans_ben = c()
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
  
  h <- hist(orfans, breaks =seq(min(insertion_index),(max(insertion_index)+1),0.02), plot = FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(2)
  par(mar=c(5.1,2.5,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab=NA, ylab=NA, main ="Unconserved", cex.axis=1.5, cex.main = 1.5,
       xlim=c(0,4), ylim=c(0,150), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,150,50), labels=c(0,NA,NA,150), cex.axis=1.5)
  text(3.1,80, paste("n =", length(orfans)), lty=1, lwd=4, cex=1.15, bty="n")
  lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  h <- hist(multiple_copies, breaks =seq(min(insertion_index),max(insertion_index)+1,0.02), plot = FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(3)
  #par(mar=c(2,1,2,1))
  par(mar=c(5.1,1,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA, ylab=NA, main ="Conserved", cex.axis = 1.5, cex.main = 1.5,
       xlim=c(0,4), ylim=c(0,150), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,150,50), labels=c(NA,NA,NA,NA), cex.axis=1.5)
  text(2.9,80, paste("n =", length(multiple_copies)), lty=1, lwd=4, cex=1.15, bty="n")
  lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  
  close.screen(all.screens = TRUE)
  
  #dev.off()
}