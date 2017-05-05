library(stringr)
library("MASS")
library("DESeq2")
library("mclust")
library(Biostrings)
fastas_dir <- "~/EnTrI/data/fasta-protein/chromosome"
genome_dir <- "~/EnTrI/data/fasta-genome/chromosome/"
plots_dir <- "~/EnTrI/data/plot-files/chromosome"
outdir <- "~/EnTrI/results/monte-carlo2/"
list_of_files <- list.files(path=plots_dir, full.names=T, recursive=FALSE)
plots = list()
plotsnw = list()
sampledplots = list(list())
sumlength = list()
numsamples = 100
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  windowsize = 100
  
  #comment
  genomefile = paste(genome_dir,locusid,'.fa', sep = '')
  g=readDNAStringSet(genomefile)[1][[1]]
  genome=strsplit(as.character(g),"")
  genomenums=(genome[[1]]=='G' | genome[[1]]=='C' | genome[[1]]=='g' | genome[[1]]=='c')+0
  genomedata = unname(tapply(genomenums, (seq_along(genomenums)-1) %/% windowsize, sum))
  
  plotsnw[[locusid]] = plotfile[,1] + plotfile[,2]
  plots[[locusid]] = unname(tapply(plotsnw[[locusid]], (seq_along(plotsnw[[locusid]])-1) %/% windowsize, sum)) # sums the number of insertions in each window
  
  # comment
  plotsdata = as.data.frame(cbind(1:length(plots[[locusid]]), as.vector(plots[[locusid]])))
  names(plotsdata) <- c("position", "num_inserts")
  mdl <- loess(log(plotsdata$num_inserts+1)~plotsdata$position, span=0.2, family = "gaussian",
               control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
  weights = exp(predict(mdl)) - 1
  # eps = 1e-10
  # weightspos = (weightspos-min(weightspos)+eps)/(max(weightspos)-min(weightspos)+eps)
  # weightspos = weightspos / sum(weightspos)
  weights = weights / median(weights)
  plots[[locusid]] = plots[[locusid]] / weights
  
  mdl <- loess(log(plots[[locusid]]+1)~genomedata, span=0.2, family = "gaussian",
              control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
  weights = exp(predict(mdl)) - 1
  # weightsgc = (weightsgc-min(weightsgc)+eps)/(max(weightsgc)-min(weightsgc)+eps)
  # weightsgc = weightsgc / sum(weightsgc)
  weights = weights / median(weights)
  plots[[locusid]] = round(plots[[locusid]] / weights)
  # weights = weightspos * weightsgc
  # weightorder = order(weights, decreasing = TRUE)
  
  for (i in 1:numsamples)
  {
    # samples = sample(plots[[locusid]], prob = weights) # The values with highest probability are more probable to appear earlier in the vector
    # sampledplots[[locusid]][[i]] = samples[weightorder]
    sampledplots[[locusid]][[i]] = sample(plots[[locusid]])
  }
  sumlength[[locusid]] = c(sum(plots[[locusid]]), length(plots[[locusid]]))
}

list_of_files <- list.files(path=fastas_dir, full.names=T, recursive=FALSE)
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/data/fasta-protein/chromosome/U00096.fasta",list_of_files) ]
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta",list_of_files) ]

for (filename in list_of_files)
{
  print(filename)
  counttable <- data.frame()
  genenames = c()
  fastafile = readLines(filename)
  iitable = list()
  for (line in fastafile)
  {
    if (startsWith(line, ">"))
    {
      matchresult = str_match(line, ">([[:graph:]]+)[[:blank:]]\\[[[:graph:]]+\\/([[:digit:]]+)\\-([[:digit:]]+)[[:print:]]+\\(([[:alpha:]]+)\\)")
      if (!(is.na(matchresult[5])))
      {
        locustag = matchresult[2]
        start = as.numeric(matchresult[3])
        end = as.numeric(matchresult[4])
        direction = matchresult[5]
        locusid = str_match(locustag, "([[:graph:]]+)\\_+[[:alnum:]]+")[2]
        if (is.na(locusid))
        {
          locusid = str_match(locustag, "([[:alpha:]]+)[[:digit:]]+")[2]
        }
        if (!(is.null(plots[[locusid]])))
        {
          len = end - start + 1
          starttrim = floor(len*5/100)
          endtrim = floor(len*20/100)
          genelength = "short"
          if (direction == "Forward")
          {
            newstart = start + starttrim
            newend = end - endtrim
            if (350 < len)
            {
              genelength = "long"
            }
          }
          else
          {
            newstart = start + endtrim
            newend = end - starttrim
            if (350 < len)
            {
              genelength = "long"
            }
          }
          newlen = newend - newstart + 1
          reads = sum(plotsnw[[locusid]][newstart:newend])
          samplestart = ceiling(newstart/windowsize)
          sampleend = floor(newend/windowsize)
          samplestartaddition = 1 - ((newstart %% windowsize) / windowsize)
          sampleendaddition = (newend %% windowsize) / windowsize
          samples = lapply(sampledplots[[locusid]], "[", (samplestart-1):(sampleend+1))
          for (i in seq(1,length(samples)))
          {
            samples[[i]][1] = round(samplestartaddition * samples[[i]][1])
            samples[[i]][length(samples[[i]])] = round(samplestartaddition * samples[[i]][length(samples[[i]])])
          }
          sampledreads = sapply(samples, sum)
          counttable = rbind(counttable,c(reads,sampledreads))
          genenames = c(genenames , locustag)
          ii = sum(plots[[locusid]][start:end]>0) / (end-start+1)
          iitable <- rbind(iitable, c(locustag, ii, genelength))
        }
      }
    }
  }
  row.names(counttable) <- genenames
  cnames <- c("Original")
  for (i in 1:numsamples)
  {
    cnames <- c(cnames, paste("Sample", as.character(i), sep=""))
  }
  colnames(counttable) <- cnames
  condition = factor(c("O",rep("S", numsamples)))
  dds <- DESeqDataSetFromMatrix(counttable, DataFrame(condition), ~ condition)
  res <- DESeq(dds)
  result<-results(res)
  
  # mod = Mclust(result$log2FoldChange, G=1:2)
  clust = kmeans(result$log2FoldChange, 2)
  clusters = clust$cluster
  if (sum(clusters==2) < sum(clusters==1))
  {
    clusters = 3 - clusters
  }
  essentiality = ifelse(clusters==1, 'Essential', 'Non-essential')
  
  # essentiality = c()
  # for (i in 1:length(result@rownames))
  # {
  #   if (result$log2FoldChange[i] > 5)
  #     #if (result$padj[i] < 0.001 & result$log2FoldChange[i] > 0)
  #     #if (mod$classification[i] == 1 & result$log2FoldChange[i] < 0)
  #   {
  #     essentiality = c(essentiality, 'Essential')
  #   }
  #   # else if (result$padj[i] < 0.001 & result$log2FoldChange[i] < 0)
  #   # {
  #   #   essentiality = c(essentiality, 'Beneficial-loss')
  #   # }
  #   else if (result$log2FoldChange[i] < -2)
  #   {
  #     essentiality = c(essentiality, 'Beneficial-loss')
  #   }
  #   else
  #   {
  #     essentiality = c(essentiality, 'Non-essential')
  #   }
  # }
  padj = result$padj / 2
  padj[result$log2FoldChange < 0] = 1 - padj[result$log2FoldChange < 0]
  
  dists = counttable
  for (i in 2:ncol(counttable))
  {
    dists[,i] = abs(dists[,i] - dists[,1])
  }
  condition = factor(c("O",rep("S", numsamples)))
  dds2 <- DESeqDataSetFromMatrix(dists, DataFrame(condition), ~ condition)
  res2 <- DESeq(dds2, fitType = 'mean')
  result2<-results(res2)
  lfc3=result2$log2FoldChange # DESeq2 LFC - distances
  lfc1=result$log2FoldChange # DESeq2 LFC 
  
  om = cbind(counttable$Original, rowMeans(counttable[,2:ncol(counttable)]))
  dists2 = cbind(om[,1],abs(om[,1]-om[,2]))
  eps = 1e-5
  lfc4 = log2((dists2[,1]+eps)/(dists2[,2]+eps)) # standard LFC - distances
  lfc2 = log2((om[,1]+eps)/(om[,2]+eps)) # standard LFC
  
  to_print = cbind(result@rownames, padj, lfc1, lfc2, lfc3, lfc4, essentiality)
  outpath = paste(outdir, locusid, ".txt", sep="")
  write.table(to_print, file=outpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  #plot(mod, what = "classification")
  #output_dir1 = '~/Desktop/tmc/'
  #outputpath = paste(output_dir1, locusid, ".txt", sep="")
  #write.table(iitable, file=outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}
