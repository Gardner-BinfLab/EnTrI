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
density = list()
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  sumlength[[locusid]] = c(sum(plots[[locusid]]>0), length(plots[[locusid]]))
  num_inserts[[locusid]] = sum(plots[[locusid]])
  density[[locusid]] = sumlength[[locusid]][[2]] / sumlength[[locusid]][[1]]
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
# hist(as.numeric(iitable$SL1344_1), breaks=0:(max(as.numeric(iitable$SL1344_1))*50+1)/50, xlim = c(0,4))
# hist(as.numeric(iitable$SL1344_3), breaks=0:(max(as.numeric(iitable$SL1344_3))*50+1)/50, xlim = c(0,4))
# hist(as.numeric(iitable$SL1344_7), breaks=0:(max(as.numeric(iitable$SL1344_7))*50+1)/50, xlim = c(0,4))
# hist(as.numeric(iitable$SL1344_5), breaks=0:(max(as.numeric(iitable$SL1344_5))*50+1)/50, xlim = c(0,4))
# hist(as.numeric(iitable$SL1344_9), breaks=0:(max(as.numeric(iitable$SL1344_9))*50+1)/50, xlim = c(0,4))

# plot(x=c(num_inserts$SL1344_1,num_inserts$SL1344_3,num_inserts$SL1344_7,num_inserts$SL1344_5,num_inserts$SL1344_9),
#      y=c(sumlength$SL1344_1[1],sumlength$SL1344_3[1],sumlength$SL1344_7[1],sumlength$SL1344_5[1],sumlength$SL1344_9[1]),
#      xlab = 'insertions', ylab='insertion sites', type='o')

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_1)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_1)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr1 = fp/(fp+tn)
essens1 = sum(res$cluster==2)
unif <- rep(0, sumlength[["SL1344_1"]][2])
ind <- sample(sumlength[["SL1344_1"]][2], sumlength[["SL1344_1"]][1], replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/(sumlength[["SL1344_1"]][1]/sumlength[["SL1344_1"]][2])
}
essthr <- max(as.numeric(iitable$SL1344_1)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif1 = fpunif/(fpunif+tnunif)
essensunif1 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_3)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_3)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr3 = fp/(fp+tn)
essens3 = sum(res$cluster==2)
unif <- rep(0, sumlength[["SL1344_3"]][2])
ind <- sample(sumlength[["SL1344_3"]][2], sumlength[["SL1344_3"]][1], replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/(sumlength[["SL1344_3"]][1]/sumlength[["SL1344_3"]][2])
}
essthr <- max(as.numeric(iitable$SL1344_3)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif3 = fpunif/(fpunif+tnunif)
essensunif3 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_7)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_7)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr7 = fp/(fp+tn)
essens7 = sum(res$cluster==2)
unif <- rep(0, sumlength[["SL1344_7"]][2])
ind <- sample(sumlength[["SL1344_7"]][2], sumlength[["SL1344_7"]][1], replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/(sumlength[["SL1344_7"]][1]/sumlength[["SL1344_7"]][2])
}
essthr <- max(as.numeric(iitable$SL1344_7)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif7 = fpunif/(fpunif+tnunif)
essensunif7 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_5)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_5)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr5 = fp/(fp+tn)
essens5 = sum(res$cluster==2)
unif <- rep(0, sumlength[["SL1344_5"]][2])
ind <- sample(sumlength[["SL1344_5"]][2], sumlength[["SL1344_5"]][1], replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/(sumlength[["SL1344_5"]][1]/sumlength[["SL1344_5"]][2])
}
essthr <- max(as.numeric(iitable$SL1344_5)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif5 = fpunif/(fpunif+tnunif)
essensunif5 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_9)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_9)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr9 = fp/(fp+tn)
essens9 = sum(res$cluster==2)
unif <- rep(0, sumlength[["SL1344_9"]][2])
ind <- sample(sumlength[["SL1344_9"]][2], sumlength[["SL1344_9"]][1], replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/(sumlength[["SL1344_9"]][1]/sumlength[["SL1344_9"]][2])
}
essthr <- max(as.numeric(iitable$SL1344_9)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif9 = fpunif/(fpunif+tnunif)
essensunif9 = sum(iisunif <= essthr)

one_rep_names <- c("SL1344_1", "SL1344_3", "SL1344_5", "SL1344_7", "SL1344_9")
one_rep_dens <- c(density$SL1344_1,density$SL1344_3[1],density$SL1344_5[1],density$SL1344_7[1],density$SL1344_9[1])
sorted_one_rep_dens <- sort(one_rep_dens)
sorted_one_rep_names <- one_rep_names[order(one_rep_dens)]

fprs = c()
fprsunif = c()
ins = c()
essens = c()
essensunif = c()
dens = c()
for (i in seq(5e3,55e4,5e3))
# for (i in seq(9,200,2))
{
  samplefrom <- min(sorted_one_rep_names[sorted_one_rep_dens<i])
  sampledplot <- plots[[samplefrom]]
  iis <- rep(0, nrow(geneinfo))
  num = round(sumlength$SL1344_10[[2]]/i)
  ind <- sample(which(sampledplot > 0), sum(sampledplot>0)-num, replace = FALSE)
  sampledplot[ind] = 0
  unif <- rep(0, length(sampledplot))
  ind <- sample(length(sampledplot), num, replace = FALSE)
  unif[ind] = 1
  iisunif <- rep(0, nrow(geneinfo))
  for (j in seq(1, length(iis)))
  {
    iis[j] <- sum(sampledplot[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
    iis[j] <- iis[j]/(num/length(sampledplot))
    iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
    iisunif[j] <- iisunif[j]/(num/length(sampledplot))
  }
  res <- dbscan(as.matrix(iis), minPts = 200, eps = 0.05)
  ess <- res$cluster[which.min(iis)]
  essthr <- max(iis[res$cluster==ess])
  fp = sum(res$cluster==ess & ecoessentiality!="1")
  fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
  tn = sum(res$cluster!=ess & ecoessentiality!="1")
  tnunif = sum(iisunif > essthr & ecoessentiality!="1")
  fpr = fp/(fp+tn)
  fprunif = fpunif/(fpunif+tnunif)
  fprs <- c(fprs, fpr)
  fprsunif <- c(fprsunif, fprunif)
  ins <- c(ins,i)
  dens <- c(dens, num)
  essens <- c(essens, sum(res$cluster==ess))
  essensunif <- c(essensunif, sum(iisunif <= essthr))
}
fprs <- append(fprs, fpr1, after = 2)
ins <- append(ins,sum(plots$SL1344_1>0), after = 2)
essens <- append(essens, essens1, after = 2)
fprsunif <- append(fprsunif, fprunif1, after = 2)
essensunif <- append(essensunif, essensunif1, after = 2)
fprs <- append(fprs, fpr3, after = 4)
ins <- append(ins, sum(plots$SL1344_3>0), after = 4)
essens <- append(essens, essens3, after = 4)
fprsunif <- append(fprsunif, fprunif3, after = 4)
essensunif <- append(essensunif, essensunif3, after = 4)
fprs <- append(fprs, fpr7, after = 46)
ins <- append(ins, sum(plots$SL1344_7>0), after = 46)
essens <- append(essens, essens7, after = 46)
fprsunif <- append(fprsunif, fprunif7, after = 46)
essensunif <- append(essensunif, essensunif7, after = 46)
fprs <- append(fprs, fpr5, after = 56)
ins <- append(ins, sum(plots$SL1344_5>0), after = 56)
essens <- append(essens, essens5, after = 56)
fprsunif <- append(fprsunif, fprunif5, after = 56)
essensunif <- append(essensunif, essensunif5, after = 56)
fprs <- append(fprs, fpr9, after = 99)
ins <- append(ins, sum(plots$SL1344_9>0), after = 99)
essens <- append(essens, essens9, after = 99)
fprsunif <- append(fprsunif, fprunif9, after = 99)
essensunif <- append(essensunif, essensunif9, after = 99)
pdf('~/EnTrI/figures/false-positive-rate_density.pdf')
par(mar=c(6.1,5.1,4.1,2.1))
plx<-predict(loess(fprs ~ ins, span = 0.2), se=T)
plot(ins,fprs, type = 'p', #ylim=c(0,0.1),
     col=ifelse(ins %in% c(sum(plots$SL1344_1>0), sum(plots$SL1344_3>0), sum(plots$SL1344_7>0), sum(plots$SL1344_5>0), sum(plots$SL1344_9>0)), "red", "dodgerblue2"),
     xlab = '# insertion sites', ylab = 'False positive rate', pch=20, cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,0.050)
     )
lines(ins,plx$fit, lwd=2)
lines(ins,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
lines(ins,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)

# plx<-predict(loess(fprsunif ~ ins, span = 0.2), se=T)
# points(ins, fprsunif)
# lines(ins,plx$fit, col=3)

plx<-predict(loess(essens ~ ins, span = 0.2), se=T)
plot(ins,essens, type = 'p', #ylim=c(0,0.1),
     col=ifelse(ins %in% c(sum(plots$SL1344_1>0), sum(plots$SL1344_3>0), sum(plots$SL1344_7>0), sum(plots$SL1344_5>0), sum(plots$SL1344_9>0)), "red", "dodgerblue2"),
     xlab = '# insertion sites', ylab = '# predicted essential genes', pch=20, cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,430)
     )
lines(ins,plx$fit, lwd=2)
lines(ins,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
lines(ins,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)

# plx<-predict(loess(essensunif ~ ins, span = 0.2), se=T)
# points(ins,essensunif)
# lines(ins,plx$fit, col=3)

dev.off()
