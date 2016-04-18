indir <- "../results/per-species-ii/"
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t")
dict = c("Citrobacter", "Escherichia coli ETEC CS17", "Enterobacter", "Escherichia coli ETEC H10407", "Escherichia coli UPEC",
         "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Salmonella enteritidis", "Salmonella typhimurium A130",
         "Salmonella typhimurium SL1344", "Salmonella typhimurium D23580", "Salmonella typhi")
names(dict) <- names
resolutions = c(11.72, 7.50, 7.74, 9.18, 11.26, 9.78, 21.76, 7.08, 19.80, 7.49, 10.54, 9.89)
names(resolutions) <- names
pdf("../results/per-species-insertion-index.pdf")
list_of_files <- list.files(path=indir, full.names=T, recursive=FALSE)
for (filename in list_of_files)
{
  iitable = read.table(filename, header = FALSE)
  m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
  temp <- split.screen(m)
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  h <- hist(iitable[,2], breaks =seq(min(iitable[,2]),max(iitable[,2])+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(1)
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "Insertion Index",
       main =dict[strsplit(basename(filename), "\\.")[[1]][1]], cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(0,NA,100), cex.axis=1.5)
  text(1.4,82.5, paste("N =", nrow(iitable)), lty=1, lwd=4, cex=1.15, bty="n",pos=4)
  text(2.4,82.5, paste("Resolution =", resolutions[strsplit(basename(filename), "\\.")[[1]][1]]), lty=1, lwd=4, cex=1.15, bty="n", pos=4)
  legend(2.4,80.5, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
         col=c("darkgoldenrod4","turquoise4", "darkmagenta"), bty="n")
  # lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  # lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
  orfans = c()
  single_copies = c()
  multiple_copies = c()
  for (i in seq(1,nrow(iitable)))
  {
    if (iitable[i,3] == 0)
    {
      orfans = c(orfans, iitable[i,2])
    }
    else if (iitable[i,3] == 1)
    {
      single_copies = c(single_copies, iitable[i,2])
    }
    else
    {
      multiple_copies = c(multiple_copies, iitable[i,2])
    }
  }
  h <- hist(orfans, breaks =seq(min(orfans),max(orfans)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(2)
  par(mar=c(5.1,2.5,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA, ylab = NA,main = "ORFan", cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(0,NA,100), cex.axis=1.5)
  
  h <- hist(single_copies, breaks =seq(min(single_copies),max(single_copies)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(3)
  par(mar=c(5.1,1,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA,ylab = NA, main ="Single copy", cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(NA,NA,NA), cex.axis=1.5)
  
  h <- hist(multiple_copies, breaks =seq(min(multiple_copies),max(multiple_copies)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  screen(4)
  par(mar=c(5.1,1,4.1,1))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA, ylab = NA, main = "Multiple copy", cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(NA,NA,NA), cex.axis=1.5)
  close.screen(all=TRUE)
  
  print(dict[strsplit(basename(filename), "\\.")[[1]][1]])
  print(length(orfans[orfans < 0.2 & orfans >= 0]))
  print(length(orfans >= 0.2))
  print(length(single_copies[single_copies<0.2 & single_copies>=0]) + length(multiple_copies[multiple_copies<0.2 & multiple_copies>=0]))
  print(length(single_copies[single_copies>=0.2]) + length(multiple_copies[multiple_copies>=0.2]))
  fresult <- fisher.test(matrix(c(length(orfans[orfans < 0.2& orfans >= 0]),length(orfans >= 0.2),length(single_copies[single_copies<0.2 & single_copies>=0]) +
                         length(multiple_copies[multiple_copies<0.2 & multiple_copies>=0]),length(single_copies[single_copies>=0.2]) +
                         length(multiple_copies[multiple_copies>=0.2])),nrow=2,byrow=TRUE))
  print(fresult$p.value)
}
dev.off()