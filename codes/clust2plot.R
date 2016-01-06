library(ggplot2)
# args <- commandArgs(trailingOnly = TRUE)
# clusters <- args[1]
clusters_path <- "../results/merge-clust-plot/final_clusters"
list_of_files <- list.files(path=clusters_path, full.names=T, recursive=FALSE)
file_II = list()
file_size = list()
for (filename in list_of_files)
{
  cluster <- as.matrix(read.table(filename))
  i_sum = 0
  l_sum = 0
  cluster_size = nrow(cluster)
  for (i in (1:cluster_size))
  {
    if (as.numeric(cluster[i,7]) >= 0)
    {
      i_sum = i_sum + as.numeric(cluster[i, 7])
      l_sum = l_sum + as.numeric(cluster[i, 8])
    }
  }
  if (l_sum > cluster_size * 60)
  {
    #file_II[basename(filename)] = i_sum / l_sum
    file_II[basename(filename)] = i_sum / cluster_size
    file_size[basename(filename)] = cluster_size
  }
}
insertion_index <- sapply(file_II, function(x){as.numeric(x[1])})
size_index <- sapply(file_size, function(x){as.numeric(x[1])})

pdf("../results/cluster-essentiality.pdf")

m <- rbind(c(0,1,0.5,1), c(0, 0.34, 0, 0.5), c(0.34, 0.67, 0, 0.5), c(0.67, 1, 0, 0.5))
temp <- split.screen(m)

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

h <- hist(insertion_index, breaks =500, plot=FALSE)
cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
screen(1)
plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "Insertion Index", main ="All clusters", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2, xlim=c(0,4), ylim=c(0,200), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,200,50), labels=c(0,NA,NA,NA,200), cex.axis=1.5)
legend(2.4,200, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4, 4),cex=1.15,
       col=c("darkgoldenrod4","turquoise4", "darkmagenta"), bty="n")
lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)

orfans = c()
orfans_non = c()
orfans_es = c()
orfans_ben = c()
single_occurrence = c()
single_non = c()
single_es = c()
single_ben = c()
multiple_copies = c()
multiple_non = c()
multiple_es = c()
multiple_ben = c()
for (item in names(size_index))
{
  if (size_index[item] <= 5)
  {
    orfans = c(orfans, insertion_index[item])
    if (insertion_index[item] < 0.2)
      orfans_es = c(orfans_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      orfans_non = c(orfans_non, insertion_index[item])
    else
      orfans_ben = c(orfans_ben, insertion_index[item])
  } 
  else if (size_index[item] <=  18)
  {
    single_occurrence = c(single_occurrence, insertion_index[item])
    if (insertion_index[item] < 0.2)
      single_es = c(single_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      single_non = c(single_non, insertion_index[item])
    else
      single_ben = c(single_ben, insertion_index[item])
  }
  else
  {
    multiple_copies = c(multiple_copies, insertion_index[item])
    if (insertion_index[item] < 0.2)
      multiple_es = c(multiple_es, insertion_index[item])
    else if (insertion_index[item] < 2)
      multiple_non = c(multiple_non, insertion_index[item])
    else
      multiple_ben = c(multiple_ben, insertion_index[item])
  }
}

h <- hist(orfans, breaks =500, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
screen(2)
par(mar=c(5.1,2.5,4.1,1))
plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab=NA, ylab=NA, main ="ORFan", cex.axis=1.5, cex.main = 1.5,
     xlim=c(0,4), ylim=c(0,150), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,150,50), labels=c(0,NA,NA,150), cex.axis=1.5)
#legend(0.13,120, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), bty="n",cex=1.5, lwd=c(4,4, 4),col=c("darkgoldenrod4","turquoise4", "darkmagenta"))
lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)

h <- hist(single_occurrence, breaks =500, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,0.18, 1.98, Inf))
screen(3)
par(mar=c(5.1,1,4.1,1))
plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab=NA, ylab=NA, main ="Single copy", cex.axis=1.5, cex.main = 1.5,
     xlim=c(0,4), ylim=c(0,150), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,150,50), labels=c(NA,NA,NA,NA), cex.axis=1.5)
#legend(0.13,120, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), bty="n",cex=1.5, lwd=c(4,4, 4),col=c("darkgoldenrod4","turquoise4", "darkmagenta"))
lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)

h <- hist(multiple_copies, breaks =500, plot = FALSE)
cuts <- cut(h$breaks, c(-Inf,0.19, 1.98, Inf))
screen(4)
#par(mar=c(2,1,2,1))
par(mar=c(5.1,1,4.1,1))
plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = NA, ylab=NA, main ="Multiple copy", cex.axis = 1.5, cex.main = 1.5,
     xlim=c(0,4), ylim=c(0,150), lty= "blank", axes=FALSE)
axis(1, at=seq(0,4,1), cex.axis=1.5)
axis(2, at=seq(0,150,50), labels=c(NA,NA,NA,NA), cex.axis=1.5)
#legend(0.13,25, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1),cex=1.5, bty="n", lwd=c(4,4, 4),col=c("darkgoldenrod4","turquoise4", "darkmagenta"))
lines(c(0.2, 0.2), c(-100,300), col = "red", lwd=3, lty = 2)
lines(c(2, 2), c(-100,300), col = "red", lwd=3, lty = 2)

close.screen(all.screens = TRUE)

dev.off()