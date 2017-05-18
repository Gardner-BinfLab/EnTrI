library(stringr)
library("dbscan")
fasta_dir <- "~/EnTrI/data/fasta-protein/chromosome/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_SL1344_FQ312003_v4.fasta"
plots_dir <- "~/EnTrI/data/plot-files/density/"
output_dir1 <- "~/EnTrI/results/density/"

list_of_files <- list.files(path=plots_dir, full.names=T, recursive=FALSE)
plots = list()
sumlength = list()
num_inserts = list()
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  sumlength[[locusid]] = c(sum(plots[[locusid]]>0), length(plots[[locusid]]))
  num_inserts[[locusid]] = sum(plots[[locusid]])
}

iitable = data.frame(locus.tag=NA, gene.length=NA ,SL1344_1=NA, SL1344_2=NA, SL1344_3=NA, SL1344_4=NA, SL1344_5=NA, SL1344_6=NA, SL1344_7=NA,
                     SL1344_8=NA, SL1344_9=NA, SL1344_10=NA)
fastafile = readLines(fasta_dir)
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
      len = end - start + 1
      starttrim = floor(len*5/100)
      endtrim = floor(len*20/100)
      genelength = "short"
      if (direction == "Forward")
      {
        newstart = start + starttrim
        newend = end - endtrim
        if (100 < len)
        {
          genelength = "long"
        }
      }
      else
      {
        newstart = start + endtrim
        newend = end - starttrim
        if (100 < len)
        {
          genelength = "long"
        }
      }
      newlen = newend - newstart + 1
      iirow <- c(locustag, genelength)
      for (locusid in names(iitable)[c(-1,-2)])
      {
        ii = (sum(plots[[locusid]][newstart:newend]>0) / newlen) / (sumlength[[locusid]][1] / sumlength[[locusid]][2])
        iirow <- c(iirow, ii)
      }
      iitable <- rbind(iitable, iirow)
    }
  }
}
iitable <- iitable[-1,]
hist(as.numeric(iitable$SL1344_1), breaks=0:(max(as.numeric(iitable$SL1344_1))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_2), breaks=0:(max(as.numeric(iitable$SL1344_2))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_3), breaks=0:(max(as.numeric(iitable$SL1344_3))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_4), breaks=0:(max(as.numeric(iitable$SL1344_4))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_5), breaks=0:(max(as.numeric(iitable$SL1344_5))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_6), breaks=0:(max(as.numeric(iitable$SL1344_6))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_7), breaks=0:(max(as.numeric(iitable$SL1344_7))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_8), breaks=0:(max(as.numeric(iitable$SL1344_8))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_9), breaks=0:(max(as.numeric(iitable$SL1344_9))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_10), breaks=0:(max(as.numeric(iitable$SL1344_10))*50+1)/50, xlim = c(0,4))
plot(x=c(num_inserts$SL1344_1,num_inserts$SL1344_2,num_inserts$SL1344_3,num_inserts$SL1344_4,num_inserts$SL1344_5,num_inserts$SL1344_6,
         num_inserts$SL1344_7,num_inserts$SL1344_8,num_inserts$SL1344_9,num_inserts$SL1344_10),
     y=c(sumlength$SL1344_1[1],sumlength$SL1344_2[1],sumlength$SL1344_3[1],sumlength$SL1344_4[1],sumlength$SL1344_5[1],sumlength$SL1344_6[1],
         sumlength$SL1344_7[1],sumlength$SL1344_8[1],sumlength$SL1344_9[1],sumlength$SL1344_10[1]), xlab = 'insertions', ylab='insertion sites'
     , type='o')

plot(x=c(num_inserts$SL1344_1,num_inserts$SL1344_3,num_inserts$SL1344_5,num_inserts$SL1344_7,num_inserts$SL1344_9),
     y=c(sumlength$SL1344_1[1],sumlength$SL1344_3[1],sumlength$SL1344_5[1],sumlength$SL1344_7[1],sumlength$SL1344_9[1]),
     xlab = 'insertions', ylab='insertion sites', type='o')

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_1)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_1)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_2)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_2)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_3)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_3)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_4)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_4)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_5)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_5)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_6)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_6)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_7)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_7)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_8)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_8)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_9)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_9)),col=res$cluster+1, pch=20)
sum(res$cluster==2)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_10)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_10)),col=res$cluster+1, pch=20)
sum(res$cluster==2)