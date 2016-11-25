library("ROCR")
library(stringr)

locus = c('BN373','ENC','ETEC','ROD','SL1344','STMMW','t','CS17','ERS227112','NCTC13441','SEN','SL3261','STM')
address = c('Klebsiella_pneumoniae_subsp_pneumoniae_Ecl8_HF536482_v1.fasta','Enterobacter_cloacae_subsp_cloacae_NCTC_9394_v1.fasta',
            'Escherichia_coli_ETEC_H10407_v2.fasta','Citrobacter_rodentium_ICC168_FN543502_v1.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhimurium_SL1344_FQ312003_v4.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_D23580_v1.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhi_Ty2_v1.fasta','CS17.fasta','Klebsiella_pneumoniae_RH201207_v0.fasta',
            'Escherichia_coli_UPEC_ST131_chromosome_v0.fasta','P125109.fasta',
            'SL3261.fasta','Salmonella_enterica_subsp_enterica_serovar_Typhimurium_A130_v0.fasta')
contingency <- array(0, dim=c(2,2,length(locus)))
dimnames(contingency)[[3]]=locus
dimnames(contingency)[[2]]=c('Essential', 'Non-essential')
dimnames(contingency)[[1]]=c('Essential', 'Non-essential')
names(dimnames(contingency))=c('Real', 'Predicted', 'Bacterium')

outdir = '../results/montecarlo-maximiseMCC/'

# contingencyvoting <- array(0, dim=c(2,2,length(locus)))
# dimnames(contingencyvoting)[[3]]=locus
# dimnames(contingencyvoting)[[2]]=c('Essential', 'Non-essential')
# dimnames(contingencyvoting)[[1]]=c('Essential', 'Non-essential')
# names(dimnames(contingencyvoting))=c('Real', 'Predicted', 'Bacterium')

for (i in seq(length(locus)))
{
  real = read.table(paste('../results/ecogenecounterparts/',locus[i],'.txt',sep = ''), as.is=TRUE, header=FALSE, sep="\t")
  names(real) <- c('gene', 'essentiality')
  
  lars = read.table(paste('../results/insertion-indices/normalised-insertion-indices-with-logodds/', locus[i], '.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  lars = lars[,c(4,3)]
  names(lars) <- c('logodds', 'essentiality')
  lars$logodds = -lars$logodds
  lars$essentiality[lars$essentiality == 'non-essential'] = 0
  lars$essentiality[lars$essentiality == 'beneficial-loss'] = 0
  lars$essentiality[lars$essentiality == 'essential'] = 1
  
  montecarlo = read.table(paste('../results/monte-carlo/',locus[i],'.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  montecarlo = montecarlo[,c(1,3,4)]
  names(montecarlo) <- c('gene', 'logfoldchange', 'essentiality')
  # montecarlo$logfoldchange = -montecarlo$logfoldchange
  montecarlo$essentiality[montecarlo$essentiality == 'Non-essential'] = 0
  montecarlo$essentiality[montecarlo$essentiality == 'Beneficial-loss'] = 0
  montecarlo$essentiality[montecarlo$essentiality == 'Essential'] = 1
  
  pdf(paste('../figures/essential-call-comparison-', locus[i], '.pdf', sep=''))
  
  
  predlars <- prediction(lars$logodds, real$essentiality)
  perflars <- performance(predlars,"tpr","fpr")
  auclars <- performance(predlars,measure = "auc")@y.values[[1]]
  plot(perflars,col="blue",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, main=locus[i])
  
  perflarsmcc <- performance(predlars,'mat')
  maxmcc = max(perflarsmcc@y.values[[1]][!is.na(perflarsmcc@y.values[[1]])])
  cutoff = perflarsmcc@x.values[[1]][perflarsmcc@y.values[[1]]==maxmcc & !is.na(perflarsmcc@y.values[[1]])]
  cutofflars = cutoff
  
  predmontecarlo <- prediction(montecarlo$logfoldchange, real$essentiality)
  perfmontecarlo <- performance(predmontecarlo,"tpr","fpr")
  aucmontecarlo <- performance(predmontecarlo,measure = "auc")@y.values[[1]]
  plot(perfmontecarlo,col="limegreen",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfmontecarlomcc <- performance(predmontecarlo,'mat')
  maxmcc = max(perfmontecarlomcc@y.values[[1]][!is.na(perfmontecarlomcc@y.values[[1]])])
  cutoff = perfmontecarlomcc@x.values[[1]][perfmontecarlomcc@y.values[[1]]==maxmcc & !is.na(perfmontecarlomcc@y.values[[1]])]
  cutoffmontecarlo = cutoff
  
  montecarlop = read.table(paste('../results/monte-carlo/',locus[i],'.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  montecarlop = montecarlop[,c(2,4)]
  names(montecarlop) <- c('pval', 'essentiality')
  montecarlop$pval = -montecarlop$pval
  montecarlop$essentiality[montecarlop$essentiality == 'Non-essential'] = 0
  montecarlop$essentiality[montecarlop$essentiality == 'Beneficial-loss'] = 0
  montecarlop$essentiality[montecarlop$essentiality == 'Essential'] = 1
  
  predmontecarlop <- prediction(montecarlop$pval, real$essentiality)
  perfmontecarlop <- performance(predmontecarlop,"tpr","fpr")
  aucmontecarlop <- performance(predmontecarlop,measure = "auc")@y.values[[1]]
  plot(perfmontecarlop,col="red",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfmontecarlopmcc <- performance(predmontecarlop,'mat')
  maxmcc = max(perfmontecarlopmcc@y.values[[1]][!is.na(perfmontecarlopmcc@y.values[[1]])])
  cutoff = perfmontecarlopmcc@x.values[[1]][perfmontecarlopmcc@y.values[[1]]==maxmcc & !is.na(perfmontecarlopmcc@y.values[[1]])]
  cutoffmontecarlop = cutoff
  
  montecarlod = read.table(paste('../results/monte-carlo/',locus[i],'.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  montecarlod = montecarlod[,c(5,4)]
  names(montecarlod) <- c('distlogFC', 'essentiality')
  montecarlod$distlogFC = montecarlod$distlogFC
  montecarlod$essentiality[montecarlod$essentiality == 'Non-essential'] = 0
  montecarlod$essentiality[montecarlod$essentiality == 'Beneficial-loss'] = 0
  montecarlod$essentiality[montecarlod$essentiality == 'Essential'] = 1
  
  predmontecarlod <- prediction(montecarlod$distlogFC, real$essentiality)
  perfmontecarlod <- performance(predmontecarlod,"tpr","fpr")
  aucmontecarlod <- performance(predmontecarlod,measure = "auc")@y.values[[1]]
  plot(perfmontecarlod,col="cyan",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfmontecarlodmcc <- performance(predmontecarlod,'mat')
  maxmcc = max(perfmontecarlodmcc@y.values[[1]][!is.na(perfmontecarlodmcc@y.values[[1]])])
  cutoff = perfmontecarlodmcc@x.values[[1]][perfmontecarlodmcc@y.values[[1]]==maxmcc & !is.na(perfmontecarlodmcc@y.values[[1]])]
  cutoffmontecarlod = cutoff
  
  # labels <- c(paste("Lars' method, AUC = ", format(round(auclars, 4), nsmall = 4)), paste("Monte Carlo, AUC = ", format(round(aucmontecarlo, 4), nsmall = 4)))
  # legend("bottomright", inset=.05, labels, lwd=2, col=c("blue","limegreen"))
  
  #dev.off()
  # Outlier:
  # >ETEC_0152 [FN649414/179928-180272 (Forward)] [conserved hypothetical protein]
  # MSDDVALPLEFTDAAANKVKSLIADEDNPNLKLRVYITGGGCSGFQYGFTFDDQVNEGDM
  # TIEKQGVGLVVDPMSLQYLVGGSVDYTEGLEGSRFIVTNPNAKSTCGCGSSFSI
  
  fastas_dir <- paste("~/EnTrI/data/fasta-protein/chromosome",address[i],sep='/')
  plots_dir <- paste("~/EnTrI/data/plot-files/chromosome/",locus[i],".plot",sep="")
  plots = list()
  sumlength = list()
  plotfile = as.matrix(read.table(plots_dir, as.is=TRUE))
  locusid = strsplit(basename(plots_dir),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  plotsdata = as.data.frame(cbind(1:length(plots[[locusid]]), as.vector(plots[[locusid]])))
  names(plotsdata) <- c("position", "num_inserts")
  mdl <- loess(num_inserts~position, plotsdata, span=0.2, family = "gaussian",
               control=loess.control(statistics=c("approximate"),trace.hat=c("approximate")))
  weights = predict(mdl)
  weights = weights / median(weights)
  plots[[locusid]] = round(plots[[locusid]] / weights)
  sumlength[[locusid]] = c(sum(plots[[locusid]]), length(plots[[locusid]]))
  avg = sumlength[[locusid]][1] / sumlength[[locusid]][2]
  
  #counttable <- data.frame()
  consecutivezeros = c()
  meandist = c()
  genenames = c()
  insites = c()
  fastafile = readLines(fastas_dir)
  iitable = list()
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
          consecutives = rle(plots[[locusid]][start:end])
          meanconsecutives = mean(consecutives$lengths[consecutives$values==0])/len
          maxconsecutivezeros = max(consecutives$lengths[consecutives$values==0])/len
          meandist = c(meandist, meanconsecutives)
          consecutivezeros = c(consecutivezeros, maxconsecutivezeros)
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
          reads = sum(plots[[locusid]][newstart:newend]>0) / newlen
          #samples = newlen * avg
          #sampledreads = sapply(samples, sum)
          insites = c(insites, reads)
        }
      }
    }
  }
  
  # predmean <- prediction(1/insites, real$essentiality)
  # perfmean <- performance(predmean,"tpr","fpr")
  # aucmean <- performance(predmean,measure = "auc")@y.values[[1]]
  # plot(perfmean,col="cyan",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  predconz <- prediction(consecutivezeros, real$essentiality)
  perfconz <- performance(predconz,"tpr","fpr")
  aucconz <- performance(predconz,measure = "auc")@y.values[[1]]
  plot(perfconz,col="orange",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfconzmcc <- performance(predconz,'mat')
  maxmcc = max(perfconzmcc@y.values[[1]][!is.na(perfconzmcc@y.values[[1]])])
  cutoff = perfconzmcc@x.values[[1]][perfconzmcc@y.values[[1]]==maxmcc & !is.na(perfconzmcc@y.values[[1]])]
  cutoffconz = cutoff
  
  predmeandist <- prediction(meandist, real$essentiality)
  perfmeandist <- performance(predmeandist,"tpr","fpr")
  aucmeandist <- performance(predmeandist,measure = "auc")@y.values[[1]]
  plot(perfmeandist,col="purple",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfmeandistmcc <- performance(predmeandist,'mat')
  maxmcc = max(perfmeandistmcc@y.values[[1]][!is.na(perfmeandistmcc@y.values[[1]])])
  cutoff = perfmeandistmcc@x.values[[1]][perfmeandistmcc@y.values[[1]]==maxmcc & !is.na(perfmeandistmcc@y.values[[1]])]
  cutoffmeandist = cutoff
  
  # predall <- prediction(log2(meandist)+log2(consecutivezeros)+log2(1/insites), real$essentiality)
  # perfall <- performance(predall,"tpr","fpr")
  # aucall <- performance(predall,measure = "auc")@y.values[[1]]
  # plot(perfall,col="black",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  labels <- c(paste("BioTraDIS, AUC = ", format(round(auclars, 4), nsmall = 4))
              , paste("Monte Carlo logFC, AUC = ", format(round(aucmontecarlo, 4), nsmall = 4))
              , paste("Monte Carlo Pval, AUC = ", format(round(aucmontecarlop, 4), nsmall = 4))
              , paste("Monte Carlo Distance, AUC = ", format(round(aucmontecarlod, 4), nsmall = 4))
              #, paste("Total insertions, AUC = ", format(round(aucmean, 4), nsmall = 4))
              , paste("Largest uninterrupted fraction, AUC = ", format(round(aucconz, 4), nsmall = 4))
              , paste("Mean distance between inserts, AUC = ", format(round(aucmeandist, 4), nsmall = 4)))
  legend("bottomright", inset=.05, labels, lwd=2, col=c("blue","limegreen", "red", "cyan", "orange", "purple"))
  
  
  plot(lars$logodds,montecarlo$logfoldchange, pch=20, col=real$essentiality+1, xlab = "-BioTraDIS logodds", ylab = "Monte Carlo logfoldchange")
  labels <- c("Essential","Non-essential")
  legend("bottomright", inset=.05, labels, pch=20, col=c("red","black"))
  
  plot(lars$logodds,montecarlop$pval, pch=20, col=real$essentiality+1, xlab = "-BioTraDIS logodds", ylab = "-Monte Carlo P-value")
  labels <- c("Essential","Non-essential")
  legend("bottomright", inset=.05, labels, pch=20, col=c("red","black"))
  
  hist(montecarlo$logfoldchange, breaks=150, xlab = "Monte Carlo logFC", main = "Histogram of Monte Carlo logFC")
  hist(-montecarlop$pval, breaks=150, xlab = "Monte Carlo P-value", main = "Histogram of Monte Carlo P-value")
  hist(log2(consecutivezeros+1e-5), breaks = 150, xlab = "log2(Largest uninterrupted fraction+eps)", main = "Histogram of largest uninterrupted fraction")
  hist(log2(meandist+1e-5), breaks = 150, xlab = "log2(Mean distance between inserts+eps)", main = "Histogram of mean distance between inserts")
  
  dev.off()
  
  essentiality = ifelse(montecarlo$logfoldchange >= cutoffmontecarlo, 'essential', 'non-essential')
  to_print = cbind(montecarlo$gene, montecarlo$logfoldchange, essentiality)
  outpath = paste(outdir, locusid, ".txt", sep="")
  write.table(to_print, file=outpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  tt = length(montecarlo$essentiality[montecarlo$logfoldchange >= cutoffmontecarlo & real$essentiality == "1"])
  ft = length(montecarlo$essentiality[montecarlo$logfoldchange >= cutoffmontecarlo & real$essentiality == "0"])
  ff = length(montecarlo$essentiality[montecarlo$logfoldchange <= cutoffmontecarlo & real$essentiality == "0"])
  tf = length(montecarlo$essentiality[montecarlo$logfoldchange <= cutoffmontecarlo & real$essentiality == "1"])
  contingency[,,locus[i]]=rbind(c(tt,tf),c(ft,ff))
  
  # votes = (montecarlo$logfoldchange>=cutoffmontecarlo & consecutivezeros>=cutoffconz & meandist>=cutoffmeandist) |
  #   (montecarlo$logfoldchange<cutoffmontecarlo & consecutivezeros>=cutoffconz & meandist>=cutoffmeandist) |
  #   (montecarlo$logfoldchange>=cutoffmontecarlo & consecutivezeros<cutoffconz & meandist>=cutoffmeandist) |
  #   (montecarlo$logfoldchange>=cutoffmontecarlo & consecutivezeros>=cutoffconz & meandist<cutoffmeandist)
  # tt = length(montecarlo$essentiality[votes & real$essentiality == "1"])
  # ft = length(montecarlo$essentiality[votes & real$essentiality == "0"])
  # ff = length(montecarlo$essentiality[!votes & real$essentiality == "0"])
  # tf = length(montecarlo$essentiality[!votes & real$essentiality == "1"])
  # contingencyvoting[,,locus[i]]=rbind(c(tt,tf),c(ft,ff))
}

pdf('../figures/essentiality-call-accuracy.pdf')
fourfoldplot(contingency, conf.level=0, std='i', space=0.2)
dev.off()

# pdf('../figures/essentiality-call-voting-accuracy.pdf')
# fourfoldplot(contingencyvoting, conf.level=0, std='i', space=0.2)
# dev.off()

# library(TStools)
# data <- cbind(lars$logodds,montecarlo$logfoldchange, montecarlop$pval, consecutivezeros, meandist)
# colnames(data) <- c("BioTraDIS","montecarlo", "montecarlop", "consecutivezeros", "meandist")
# nemenyi(data,conf.int=0.95,plottype="vline")


library("rgp")
functionSet1 <- functionSet("+", "*", "-","exp")
inputVariableSet1 <- inputVariableSet("montecarlo","consesutivezeros","meandist")
constantFactorySet1 <- constantFactorySet(function() rep(rnorm(1),length(real$essentiality)))
interval1 <- as.numeric(lars$logodds)
interval2 <- as.numeric(montecarlo$logfoldchange)
interval3 <- as.numeric(consecutivezeros)
interval4 <- as.numeric(meandist)
#fitnessFunction1 <- function(f) print(f)
fitnessFunction1 <- function(f) 1-performance(prediction(f(interval2,interval3,interval4), real$essentiality),measure = "auc")@y.values[[1]]
# fitnessFunction1 <- function(f) rmse(f(interval2,interval3,interval4), as.numeric(real$essentiality))
set.seed(1)
gpResult1 <- geneticProgramming(functionSet = functionSet1, inputVariables = inputVariableSet1, constantSet = constantFactorySet1,
                                fitnessFunction = fitnessFunction1, stopCondition = makeTimeStopCondition(5 * 60))
bestSolution1 <- gpResult1$population[[which.min(gpResult1$fitnessValues)]]

predgp <- prediction(bestSolution1(interval2,interval3,interval4), real$essentiality)
perfgp <- performance(predgp,"tpr","fpr")
aucgp <- performance(predgp,measure = "auc")@y.values[[1]]
plot(perfgp,col="purple",lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7)
format(round(aucgp, 4), nsmall = 4)
