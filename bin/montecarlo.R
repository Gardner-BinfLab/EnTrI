library(stringr)
library("MASS")
library("DESeq2")
#library("mclust")
fastas_dir <- "~/EnTrI/data/fasta-protein/chromosome"
plots_dir <- "~/EnTrI/data/plot-files/chromosome"
list_of_files <- list.files(path=plots_dir, full.names=T, recursive=FALSE)
plots = list()
sampledplots = list(list())
sumlength = list()
weights = list()
numsamples = 2
locusid = c()
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  lid = strsplit(basename(filename),"\\.")[[1]][1]
  locusid = c(locusid,lid)
  plots[[lid]] = plotfile[,1] + plotfile[,2]
  # plot(plots[[lid]], pch = '.', ylim=c(0,10))
  # lines(loess.smooth(1:length(plots[[lid]]), plots[[lid]], span=0.2, family = "gaussian"), col=2, lwd=5)
  plotsdata = as.data.frame(cbind(1:length(plots[[lid]]), as.vector(plots[[lid]])))
  names(plotsdata) <- c("position", "num_inserts")
  mdl <- loess(num_inserts~position, plotsdata, span=0.2, family = "gaussian")
  # mdl <- loess.smooth(plotsdata$position, plotsdata$num_inserts, span=0.2)
  weights[[lid]] = predict(mdl)
  weights[[lid]] = weights[[lid]] / sum(weights[[lid]])
  sumlength[[lid]] = c(sum(plots[[lid]]), length(plots[[lid]]))
}
sampledplots = lapply(locusid, function(filename) lapply( 1:numsamples, function(i) sample(plots[[filename]], prob=weights[[lid]])))
names(sampledplots) <- locusid

list_of_files <- list.files(path=fastas_dir, full.names=T, recursive=FALSE)
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/data/fasta-protein/chromosome/U00096.fasta",list_of_files) ]
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta",list_of_files) ]

for (filename in list_of_files)
{
  counttable <- data.frame()
  genenames = c()
  fastafile = readLines(filename)
  for (line in fastafile)
  {
    if (startsWith(line, ">"))
    {
      matchresult = str_match(line, ">([[:graph:]]+)[[:blank:]]\\[[[:graph:]]+\\/([[:digit:]]+)\\-([[:digit:]]+)[[:blank:]]\\(([[:alpha:]]+)\\)")
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
          reads = sum(plots[[locusid]][newstart:newend])
          samples = lapply(sampledplots[[locusid]], "[", newstart:newend)
          sampledreads = sapply(samples, sum)
          counttable = rbind(counttable,c(reads,sampledreads))
          genenames = c(genenames , locustag)
          # ii = sum(plots[[locusid]][start:end]>0) / (end-start+1)
          # iitable <- rbind(iitable, c(locustag, ii, genelength))
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
  print(sum(result$padj < 0.001 & result$log2FoldChange < 0))
  # mod = Mclust(result$log2FoldChange, G=2)
  # plot(mod, what = "classification")
  #outputpath = paste(output_dir1, locusid, ".txt", sep="")
  #write.table(iitable, file=outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}