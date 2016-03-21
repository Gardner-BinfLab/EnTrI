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
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  h <- hist(iitable[,2], breaks =seq(min(iitable[,2]),max(iitable[,2])+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "Insertion Index",
       main =dict[strsplit(basename(filename), "\\.")[[1]][1]], cex.lab = 2,
       cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,100), lty= "blank", axes=FALSE)
  axis(1, at=seq(0,4,1), cex.axis=1.5)
  axis(2, at=seq(0,100,50), labels=c(0,50,100), cex.axis=1.5)
  text(2.4,87, paste("N =", nrow(iitable)), lty=1, lwd=4, cex=1.15, bty="n",pos=4)
  text(2.4,82.5, paste("Resolution =", resolutions[strsplit(basename(filename), "\\.")[[1]][1]]), lty=1, lwd=4, cex=1.15, bty="n", pos=4)
  legend(2.4,80, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
         col=c("darkgoldenrod4","turquoise4", "darkmagenta"), bty="n")
  # lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
  # lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)
}
dev.off()