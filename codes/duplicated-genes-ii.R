require(stringr)
clustpath <- "~/EnTrI/results/merge-clust-plot/final_clusters/"
copy1_ii <- c()
copy2_ii <- c()
list_of_files <- list.files(path=clustpath, pattern="*.txt", full.names=F, recursive=FALSE)
for (filename in list_of_files)
{
  clust_tbl <- read.table(paste(clustpath,'/',filename,sep = ''), sep = '\t', header = FALSE)
  for (i in seq(1,nrow(clust_tbl)-1))
  {
    if (nrow(clust_tbl) > 1)
    {
      for (j in seq(i+1, nrow(clust_tbl)))
      {
        ii1 = as.numeric(clust_tbl[i,7])
        ii2 = as.numeric(clust_tbl[j,7])
        match1 = str_match(clust_tbl[i,2], "([[:graph:]]+)\\_[[:digit:]]+")[2]
        if (is.na(match1))
        {
          match1 = str_match(clust_tbl[i,2], "([[:alpha:]]+)[[:digit:]]+")[2]
        }
        match2 = str_match(clust_tbl[j,2], "([[:graph:]]+)\\_[[:digit:]]+")[2]
        if (is.na(match2))
        {
          match2 = str_match(clust_tbl[j,2], "([[:alpha:]]+)[[:digit:]]+")[2]
        }
        if (ii1 != -1 & ii2 != -1 & match1 == match2)
        {
          maxii <- max(ii1,ii2)
          minii <- min(ii1,ii2)
          copy1_ii <- c(copy1_ii, maxii)
          copy2_ii <- c(copy2_ii, minii)
        }
      }
    }
  }
}

pdf("../results/copy1-vs-copy2_ii.pdf")
smoothScatter(copy2_ii, copy1_ii, nbin=1000, nrpoints=0, xlim=c(0,4), ylim=c(0,4))
dev.off()
