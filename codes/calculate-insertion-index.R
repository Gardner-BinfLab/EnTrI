library(stringr)
library("MASS")
fastas_dir <- "~/EnTrI/sequences/fasta-protein/chromosome"
plots_dir <- "~/EnTrI/sequences/plot-files/chromosome"
output_dir1 <- "~/EnTrI/results/insertion-indices/insfree/"
output_dir2 <- "~/EnTrI/results/insertion-indices/gamma/"
list_of_files <- list.files(path=plots_dir, full.names=T, recursive=FALSE)
plots = list()
sumlength = list()
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  sumlength[[locusid]] = c(sum(plots[[locusid]]>0), length(plots[[locusid]]))
}

list_of_files <- list.files(path=fastas_dir, full.names=T, recursive=FALSE)
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/U00096.fasta",list_of_files) ]
list_of_files <- list_of_files[ !grepl("/home/fatemeh/EnTrI/sequences/fasta-protein/chromosome/seqdb.fasta",list_of_files) ]
dir.create(output_dir1)
dir.create(output_dir2)
for (filename in list_of_files)
{
  iitable = c()
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
          ii = (sum(plots[[locusid]][start:end]>0) / len) / (sumlength[[locusid]][1] / sumlength[[locusid]][2])
          # ii = sum(plots[[locusid]][start:end]>0) / (end-start+1)
          essentiality = "not"
          if (direction == "Forward")
          {
            if (100 < len & sum(plots[[locusid]][(start+starttrim):(end-endtrim)]) == 0)
            {
              essentiality = "essential"
            }
            iitable <- rbind(iitable, c(locustag, ii, essentiality))
          }
          else
          {
            if (100 < len & sum(plots[[locusid]][(start+endtrim):(end-starttrim)]) == 0)
            {
              essentiality = "essential"
            }
            iitable <- rbind(iitable, c(locustag, ii, essentiality))
          }
        }
      }
    }
  }
  outputpath = paste(output_dir1, locusid, ".txt", sep="")
  write.table(iitable, file=outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}

names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names

pdf("~/EnTrI/results/insertion-indices/results.pdf")
list_of_files <- list.files(path=output_dir1, full.names=T, recursive=FALSE)
for (filename in list_of_files)
{
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  iifile = as.matrix(read.table(filename, as.is = TRUE))
  
  ii= as.numeric(iifile[,2])
  nG = length(ii)
  
  #identify second maxima
  h <- hist(ii, breaks=100,plot=FALSE)
  maxindex <- which.max(h$density[3:length(h$density)])
  maxval <- h$mids[maxindex+2]
  
  #find inter-mode minima with loess
  r <- floor(maxval *1000)
  I = ii < r / 1000
  h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
  lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  # plot(h1$density, main="Density",ylim = c(0,5))
  # lines(predict(lo),col='red',lwd=2)
  m1 = h1$mids[which.min(predict(lo))]
  m2 = h$mids[max(which(h$counts>5))]
  I1 = ((ii < m1)&(ii > 0))
  I2 = ((ii >= m1)&(ii < m2))
  # I3 = (ii >= m2)
  
  f1 = (sum(I1) + sum(ii == 0))/nG
  f2 = (sum(I2))/nG
  # f3 = (sum(I3))/nG
  
  d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
  d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves
  # d3 = fitdistr(ii[I3], "gamma")

  #plots
  hist(ii,breaks=200, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  lines(0:2000/500, f1*dgamma(0:2000/500, 1, d1$estimate[2])) # was [2]
  lines(0:2000/500, f2*dgamma(0:2000/500, d2$estimate[1], d2$estimate[2]))
  # print changepoint
  
  #calculate log-odds ratios to choose thresholds
  lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
  upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
  
  essen <- lower/10000
  ambig <- upper/10000
  # noness <- m2
  noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  
  lines(c(essen, essen), c(0,20), col="red")
  lines(c(ambig, ambig), c(0,20), col="red")
  lines(c(noness, noness), c(0,20), col="red")
  
  mtext(paste(essen, ":", "Essential changepoint"), side=3, adj=1, padj=2)
  mtext(paste(ambig, ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
  mtext(paste(noness, ":", "Non-essential changepoint"), side=3, adj=1, padj=5.5)
  ###https://cran.r-project.org/web/views/Cluster.html
  for (i in (1:length(iifile[,1])))
  {
    if (iifile[i,2] < essen)
    {
      iifile[i,3] = "essential"
    }
    else if(iifile[i,2] < noness)
    {
      iifile[i,3] = "non-essential"
    }
    else
    {
      iifile[i,3] = "beneficial-loss"
    }
  }
  outputpath = filename
  outputpath = gsub("insfree", "gamma", filename)
  write.table(iifile, file=outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}
dev.off()