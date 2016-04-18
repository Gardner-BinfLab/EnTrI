ratios_conserved <- c(8*100/(320+8), 31*100/(135+31), 36*100/(159+36), 11*100/(240+11), 14*100/(186+14), 24*100/(582+21), 100*100/(488+100), 48*100/(600+48), 40*100/(791+40), 40*100/(1228+40), 30*100/(711+30), 41*100/(604+41))
ratios_unconserved <- c(385*100/(3583+385), 363*100/(3772+363), 376*100/(3968+376), 427*100/(3767+427), 361*100/(3934+361), 354*100/(3800+354), 395*100/(3677+395), 374*100/(3728+374), 371*100/(3783+371), 376*100/(4018+376), 305*100/(3949+305), 267*100/(2811+267))
ratios <- rbind(ratios_conserved,ratios_unconserved)
p_values <- c(7.323e-07, 8.904e-05, 3.020307e-05, 0.001477236, 0.6000692, 3.70441e-05, 4.325686e-07, 0.1807112, 2.494271e-05, 2.776819e-12, 0.001383878, 0.05882195)
names <- c("Salmonella typhi", "Salmonella enteritidis", "Salmonella typhimurium SL1344", "Salmonella typhimurium A130",
           "Salmonella typhimurium D23580", "Escherichia coli UPEC", "Escherichia coli ETEC CS17", "Escherichia coli ETEC H10407",
           "Citrobacter", "Klebsiella pneumoniae RH201207", "Klebsiella pneumoniae Ecl8", "Enterobacter")
names(p_values) <- names
colnames(ratios) <- names
mar.default <- c(5,4,2,2) + 0.1
par(mar = mar.default + c(13, 1, 0, 0))
barplot(ratios, beside=TRUE,ylim=c(0,25),col=c("limegreen","gray17"), las=2, ylab="% Essential", cex.lab = 1.5, axes=FALSE)
axis(2, at=c(0,25), cex.axis=1.5)
text(1.8, max(ratios[,1])+1, "***", cex=2)
text(4.7, max(ratios[,2])+1, "***", cex=2)
text(7.75, max(ratios[,3])+1, "***", cex=2)
text(10.7, max(ratios[,4])+1, "*", cex=2)
text(13.75, max(ratios[,5])+1.5, "NS", cex=1)
text(16.7, max(ratios[,6])+1, "***", cex=2)
text(19.7, max(ratios[,7])+1, "***", cex=2)
text(22.7, max(ratios[,8])+1.5, "NS", cex=1)
text(25.7, max(ratios[,9])+1, "***", cex=2)
text(28.75, max(ratios[,10])+1, "***", cex=2)
text(31.75, max(ratios[,11])+1, "**", cex=2)
text(34.7, max(ratios[,12])+1.5, "NS", cex=1)
legend(25,25, legend=c("Unconserved","Conserved"), lty=c(1,1), lwd=c(4,4),cex=1.15, col=c("limegreen","gray17"), bty="n")