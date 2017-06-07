library(stringr)
library("dbscan")
fasta_dir <- "~/EnTrI/data/fasta-protein/chromosome/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_SL1344_FQ312003_v4.fasta"
plots_dir <- "~/EnTrI/data/plot-files/density/"
output_dir1 <- "~/EnTrI/results/density/"
ecogene <- "~/EnTrI/results/ecogenecounterparts/SL1344.txt"

ecoessentiality <- read.table(ecogene, row.names = 1, col.names = c('name','essentiality'))

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
geneinfo <- matrix(ncol=3)
colnames(geneinfo) <- c('locustag', 'start', 'end')
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
      geneinfo <- rbind(geneinfo, c(locustag, newstart, newend))
    }
  }
}
iitable <- iitable[-1,]
geneinfo <- geneinfo[-1,]
hist(as.numeric(iitable$SL1344_1), breaks=0:(max(as.numeric(iitable$SL1344_1))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_3), breaks=0:(max(as.numeric(iitable$SL1344_3))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_7), breaks=0:(max(as.numeric(iitable$SL1344_7))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_5), breaks=0:(max(as.numeric(iitable$SL1344_5))*50+1)/50, xlim = c(0,4))
hist(as.numeric(iitable$SL1344_9), breaks=0:(max(as.numeric(iitable$SL1344_9))*50+1)/50, xlim = c(0,4))
# plot(x=c(num_inserts$SL1344_1,num_inserts$SL1344_2,num_inserts$SL1344_3,num_inserts$SL1344_4,num_inserts$SL1344_5,num_inserts$SL1344_6,
#          num_inserts$SL1344_7,num_inserts$SL1344_8,num_inserts$SL1344_9,num_inserts$SL1344_10),
#      y=c(sumlength$SL1344_1[1],sumlength$SL1344_2[1],sumlength$SL1344_3[1],sumlength$SL1344_4[1],sumlength$SL1344_5[1],sumlength$SL1344_6[1],
#          sumlength$SL1344_7[1],sumlength$SL1344_8[1],sumlength$SL1344_9[1],sumlength$SL1344_10[1]), xlab = 'insertions', ylab='insertion sites'
#      , type='o')



plot(x=c(num_inserts$SL1344_1,num_inserts$SL1344_3,num_inserts$SL1344_7,num_inserts$SL1344_5,num_inserts$SL1344_9),
     y=c(sumlength$SL1344_1[1],sumlength$SL1344_3[1],sumlength$SL1344_7[1],sumlength$SL1344_5[1],sumlength$SL1344_9[1]),
     xlab = 'insertions', ylab='insertion sites', type='o')

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_1)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_1)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr1 = fp/(fp+tn)
print(fpr1)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_3)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_3)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr3 = fp/(fp+tn)
print(fpr3)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_7)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_7)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr7 = fp/(fp+tn)
print(fpr7)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_5)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_5)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr5 = fp/(fp+tn)
print(fpr5)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_9)), minPts = 200, eps = 0.05)
plot(as.matrix(as.numeric(iitable$SL1344_9)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr9 = fp/(fp+tn)
print(fpr9)

one_rep_names <- c("SL1344_1", "SL1344_3", "SL1344_5", "SL1344_7", "SL1344_9")
one_rep_inserts <- c(num_inserts$SL1344_1,num_inserts$SL1344_3,num_inserts$SL1344_5,num_inserts$SL1344_7,num_inserts$SL1344_9)
sorted_one_rep_inserts <- sort(one_rep_inserts)
sorted_one_rep_names <- one_rep_names[order(one_rep_inserts)]

fprs = c()
ins = c()
for (i in seq(4e6,138e5,1e5))
{
  samplefrom <- min(sorted_one_rep_names[sorted_one_rep_inserts>i])
  sampledplot <- plots[[samplefrom]]
  iis <- rep(0, nrow(geneinfo))
  while (sum(sampledplot) > i)
  {  
    ind <- sample(which(sampledplot > 0), (sum(sampledplot)-i), replace = TRUE)
    tbl <- table(ind)
    for (item in names(tbl))
    {
      sampledplot[as.numeric(item)] = max((sampledplot[as.numeric(item)] - tbl[item]), 0)
    }
  }
  for (j in seq(1, length(iis)))
  {
    iis[j] <- sum(sampledplot[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])])/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2]))
    iis[j] <- iis[j]/(i/length(sampledplot))
  }
  res <- dbscan(as.matrix(iis), minPts = 200, eps = 0.05)
  ess <- res$cluster[which.min(iis)]
  fp = sum(res$cluster==ess & ecoessentiality!="1")
  tn = sum(res$cluster!=ess & ecoessentiality!="1")
  fpr = fp/(fp+tn)
  fprs <- c(fprs, fpr)
  ins <- c(ins,i)
  print(i)
  print(fpr)
}
fprs <- append(fprs, fpr3, after = 1)
ins <- append(ins,sum(plots$SL1344_3), after = 1)
fprs <- append(fprs, fpr1, after = 3)
ins <- append(ins, sum(plots$SL1344_1), after = 3)
fprs <- append(fprs, fpr5, after = 8)
ins <- append(ins, sum(plots$SL1344_5), after = 8)
fprs <- append(fprs, fpr7, after = 21)
ins <- append(ins, sum(plots$SL1344_7), after = 21)
fprs <- append(fprs, fpr9, after = 103)
ins <- append(ins, sum(plots$SL1344_9), after = 103)
pdf('~/EnTrI/figures/false-positive-rate_density.pdf')
plot(ins,fprs, type = 'o', #ylim=c(0,0.1),
     col=ifelse(ins %in% c(sum(plots$SL1344_3), sum(plots$SL1344_1), sum(plots$SL1344_5), sum(plots$SL1344_7), sum(plots$SL1344_9)), "red", "black"),
     xlab = 'number of insertions', ylab = 'False positive rate')
dev.off()