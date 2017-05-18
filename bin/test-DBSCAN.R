library("dbscan")

genome_l = 5000001
gene_n = 5000
no_ins = 230
# insertion_n = round(genome_l/50)
pdf('~/EnTrI/figures/simulation.pdf')
for (insertion_n in round(genome_l/seq(10,100,10)))
{
  coordinates <- sort(sample(seq(1,genome_l/3), gene_n*2))
  starts <- coordinates[1:2%%2==1]*3
  ends <- (coordinates[1:2%%2==0]*3 + 51) %% genome_l
  if (ends[length(ends)] == 0)
  {
    ends[length(ends)] = genome_l
  }
  essentiality <- sample(c(0,1), gene_n, replace = TRUE, prob = c((gene_n-no_ins)/gene_n, no_ins/gene_n))
  essentials_l = sum(ends[essentiality==1]-starts[essentiality==1]+1)
  insertions <- table(sample(seq(1,genome_l-essentials_l), insertion_n, replace = TRUE))
  plot <- c(rep(0, genome_l-essentials_l))
  plot[as.numeric(names(insertions))] = insertions[names(insertions)]
  essential_starts <- starts[essentiality==1]
  essential_ends <- ends[essentiality==1]
  for (i in seq(1,length(essential_starts)))
  {
    leni = essential_ends[i]-essential_starts[i]+ 1
    toappend <- rep(0,leni)
    start = essential_starts[i] - 1
    plot <- append(plot, toappend, start)
  }
  insertion_indices <- rep(0,gene_n)
  for (i in seq(1,length(starts)))
  {
    if (starts[i] < ends[i])
    {
      insertion_indices[i] = (sum(plot[starts[i]:ends[i]]) / (ends[i] - starts[i] + 1)) / (insertion_n/genome_l)
    }else
    {
      insertion_indices[i] = ((sum(plot[starts[i]:genome_l])+sum(plot[1:ends[i]])) / (ends[i] - starts[i] + 1)) / (insertion_n/genome_l)
    }
  }
  res <- dbscan(as.matrix(insertion_indices), minPts = 200, eps = 0.05)
  plot(insertion_indices,col=res$cluster+1, main=paste('Simulated data -', toString(insertion_n), 'insertions'), xlab='Gene', ylab='Insertion index')
  ess <- res$cluster[which.min(insertion_indices)]
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  nes <- getmode(res$cluster)
  belthr <- max(insertion_indices[res$cluster==nes])
  essthr <- max(insertion_indices[res$cluster==ess])
  nesthr <- min(insertion_indices[res$cluster==nes])
  hist(insertion_indices,breaks=0:(max(insertion_indices)*100+1)/100, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=paste('Simulated data -', toString(insertion_n), 'insertions'))
  text(2,4,paste(toString(sum(res$cluster == ess)), "essential genes"))
  lines(c(essthr, essthr), c(0,20), col="red")
}
dev.off()