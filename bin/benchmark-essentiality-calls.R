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
colors=c('blue', 'limegreen', 'red', 'cyan', 'black', 'orange', 'purple', 'gray')
for (i in seq(length(locus)))
{
  real = read.table(paste('../results/ecogenecounterparts/',locus[i],'.txt',sep = ''), as.is=TRUE, header=FALSE, sep="\t")
  names(real) <- c('gene', 'essentiality')
  
  biotradis = read.table(paste('../results/insertion-indices/normalised-insertion-indices-with-logodds/', locus[i], '.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  biotradis = biotradis[,4]
  biotradis = -biotradis
  
  montecarlo = read.table(paste('../results/monte-carlo/',locus[i],'.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  montecarlo = montecarlo[,c(2,3,4,5,6)]
  montecarlo[,1] = -montecarlo[,1]
  montecarlo[,3] = -montecarlo[,3]
  montecarlo[,5] = -montecarlo[,5]
  names(montecarlo) <- c('-DESeqPval','DESeqLFC','-LFC','DESeqLFCDist', '-LFCDist')
  # montecarlo$logfoldchange = -montecarlo$logfoldchange
  
  pdf(paste('../figures/essential-call-comparison-', locus[i], '.pdf', sep=''))
  
  predbiotradis <- prediction(biotradis, real$essentiality)
  perfbiotradis <- performance(predbiotradis,"tpr","fpr")
  aucbiotradis <- performance(predbiotradis,measure = "auc")@y.values[[1]]
  plot(perfbiotradis,col=colors[1],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, main=locus[i])
  
  perfbiotradismcc <- performance(predbiotradis,'mat')
  maxmcc = max(perfbiotradismcc@y.values[[1]][!is.na(perfbiotradismcc@y.values[[1]])])
  cutoff = perfbiotradismcc@x.values[[1]][perfbiotradismcc@y.values[[1]]==maxmcc & !is.na(perfbiotradismcc@y.values[[1]])]
  cutoffbiotradis = cutoff
  
  aucmontecarlo = c()
  cutoffmontecarlo = c()
  for (j in seq(5))
  {
    predmontecarlo <- prediction(montecarlo[,j], real$essentiality)
    perfmontecarlo <- performance(predmontecarlo,"tpr","fpr")
    aucmontecarlo <- c(aucmontecarlo, performance(predmontecarlo,measure = "auc")@y.values[[1]])
    plot(perfmontecarlo,col=colors[j+1],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
    
    perfmontecarlomcc <- performance(predmontecarlo,'mat')
    maxmcc = max(perfmontecarlomcc@y.values[[1]][!is.na(perfmontecarlomcc@y.values[[1]])])
    cutoff = perfmontecarlomcc@x.values[[1]][perfmontecarlomcc@y.values[[1]]==maxmcc & !is.na(perfmontecarlomcc@y.values[[1]])]
    cutoffmontecarlo = c(cutoffmontecarlo,cutoff)
  }
  
  fastas_dir <- paste("~/EnTrI/data/fasta-protein/chromosome",address[i],sep='/')
  plots_dir <- paste("~/EnTrI/data/plot-files/chromosome/",locus[i],".plot",sep="")
  plots = list()
  sumlength = list()
  plotfile = as.matrix(read.table(plots_dir, as.is=TRUE))
  locusid = strsplit(basename(plots_dir),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  
  
  #counttable <- data.frame()
  consecutivezeros = c()
  meandist = c()
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
          meanconsecutives = mean(consecutives$lengths[consecutives$values==0])
          maxconsecutivezeros = max(consecutives$lengths[consecutives$values==0])
          meandist = c(meandist, meanconsecutives)
          consecutivezeros = c(consecutivezeros, maxconsecutivezeros)
        }
      }
    }
  }
  predconz <- prediction(consecutivezeros, real$essentiality)
  perfconz <- performance(predconz,"tpr","fpr")
  aucconz <- performance(predconz,measure = "auc")@y.values[[1]]
  plot(perfconz,col=colors[7],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfconzmcc <- performance(predconz,'mat')
  maxmcc = max(perfconzmcc@y.values[[1]][!is.na(perfconzmcc@y.values[[1]])])
  cutoff = perfconzmcc@x.values[[1]][perfconzmcc@y.values[[1]]==maxmcc & !is.na(perfconzmcc@y.values[[1]])]
  cutoffconz = cutoff
  
  predmeandist <- prediction(meandist, real$essentiality)
  perfmeandist <- performance(predmeandist,"tpr","fpr")
  aucmeandist <- performance(predmeandist,measure = "auc")@y.values[[1]]
  plot(perfmeandist,col=colors[8],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  
  perfmeandistmcc <- performance(predmeandist,'mat')
  maxmcc = max(perfmeandistmcc@y.values[[1]][!is.na(perfmeandistmcc@y.values[[1]])])
  cutoff = perfmeandistmcc@x.values[[1]][perfmeandistmcc@y.values[[1]]==maxmcc & !is.na(perfmeandistmcc@y.values[[1]])]
  cutoffmeandist = cutoff
  
  labels <- c(paste("BioTraDIS, AUC = ", format(round(aucbiotradis, 4), nsmall = 4))
              , paste("Monte Carlo Pval, AUC = ", format(round(aucmontecarlo[1], 4), nsmall = 4))
              , paste("Monte Carlo DESeq LFC, AUC = ", format(round(aucmontecarlo[2], 4), nsmall = 4))
              , paste("Monte Carlo LFC, AUC = ", format(round(aucmontecarlo[3], 4), nsmall = 4))
              , paste("Monte Carlo DESeq LFC Distances, AUC = ", format(round(aucmontecarlo[4], 4), nsmall = 4))
              , paste("Monte Carlo LFC Distances, AUC = ", format(round(aucmontecarlo[5], 4), nsmall = 4))
              , paste("Largest uninterrupted fraction, AUC = ", format(round(aucconz, 4), nsmall = 4))
              , paste("Mean distance between inserts, AUC = ", format(round(aucmeandist, 4), nsmall = 4)))
  legend("bottomright", inset=.05, labels, lwd=2, col=colors)
  
  plot(-biotradis,montecarlo$DESeqLFC, pch=20, col=real$essentiality+1, xlab = "BioTraDIS logodds", ylab = "Monte Carlo logfoldchange")
  labels <- c("Essential","Non-essential")
  legend("topright", inset=.05, labels, pch=20, col=c("red","black"))
  
  plot(-biotradis,log2(-montecarlo$`-DESeqPval`), pch=20, col=real$essentiality+1, xlab = "BioTraDIS logodds", ylab = "log2(Monte Carlo P-value)")
  labels <- c("Essential","Non-essential")
  legend("bottomright", inset=.05, labels, pch=20, col=c("red","black"))
  
  hist(montecarlo$DESeqLFC, breaks=150, xlab = "Monte Carlo logFC", main = "Histogram of Monte Carlo logFC")
  hist(log2(-montecarlo$`-DESeqPval`), breaks=150, xlab = "log2(Monte Carlo P-value)", main = "Histogram of Monte Carlo P-value")
  hist(log2(consecutivezeros+1e-5), breaks = 150, xlab = "log2(Largest uninterrupted fraction+1e-5)", main = "Histogram of largest uninterrupted fraction")
  hist(log2(meandist+1e-5), breaks = 150, xlab = "log2(Mean distance between inserts+1e-5)", main = "Histogram of mean distance between inserts")
  
  dev.off()
  
  essentiality = ifelse(montecarlo$DESeqLFC >= cutoffmontecarlo[2], 'essential', 'non-essential')
  to_print = cbind(real$gene, montecarlo$DESeqLFC, essentiality)
  outpath = paste(outdir, locusid, ".txt", sep="")
  write.table(to_print, file=outpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  tt = length(real$gene[montecarlo$DESeqLFC >= cutoffmontecarlo[2] & real$essentiality == "1"])
  ft = length(real$gene[montecarlo$DESeqLFC >= cutoffmontecarlo[2] & real$essentiality == "0"])
  ff = length(real$gene[montecarlo$DESeqLFC <= cutoffmontecarlo[2] & real$essentiality == "0"])
  tf = length(real$gene[montecarlo$DESeqLFC <= cutoffmontecarlo[2] & real$essentiality == "1"])
  contingency[,,locus[i]]=rbind(c(tt,tf),c(ft,ff))
  # if (locus[i] == 'CS17')
  # {
  #   cs17fps = real$gene[montecarlo$DESeqLFC >= cutoffmontecarlo[2] & real$essentiality == "0"]
  #   write.table(cs17fps, file='../results/cs17fps.txt
  #               ', quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  # }
}

pdf('../figures/essentiality-call-accuracy.pdf')
fourfoldplot(contingency, conf.level=0, std='i', space=0.2)
dev.off()