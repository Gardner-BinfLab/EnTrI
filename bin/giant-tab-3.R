library(UpSetR)
tabpath = '~/Dropbox/EnTrI/results/giant-tab/giant-tab-all.tsv'
table = read.csv(tabpath, sep='\t', row.names = 1)[,75:98]

pdf('~/Dropbox/EnTrI/figures/taxonomy-upsetr.pdf')
for (thresh in c(80,90,100))
{
  ########## with symbionts
  name = paste('Enterobacteriaceae',thresh,sep='_')
  ent = rownames(table)[table[name]==1]
  name = paste('Endosymbiont',thresh,sep='_')
  end = rownames(table)[table[name]==1]
  name = paste('Gammaproteobacter',thresh,sep='_')
  gam = rownames(table)[table[name]==1]
  name = paste('Proteobacter',thresh,sep='_')
  prot = rownames(table)[table[name]==1]
  name = paste('Bacteria',thresh,sep='_')
  bact = rownames(table)[table[name]==1]
  inpt = list(Enterobacteriaceae = ent, Endosymbiont = end, Gammaproteobacteria = gam, Proteobacteria = prot, Bacteria = bact)
  upset(fromList(inpt), order.by = "freq", mainbar.y.label = paste("Core threshold = ",thresh,"%",sep=''))
  ########## without symbionts
  name = paste('Enterobacteriaceae',thresh,sep='_')
  ent = rownames(table)[table[name]==1]
  name = paste('Gammaproteobacter.no.symb',thresh,sep='_')
  gam = rownames(table)[table[name]==1]
  name = paste('Proteobacter.no.symb',thresh,sep='_')
  prot = rownames(table)[table[name]==1]
  name = paste('Bacteria.no.symb',thresh,sep='_')
  bact = rownames(table)[table[name]==1]
  inpt = list(Enterobacteriaceae = ent, Gammaproteobacteria = gam, Proteobacteria = prot, Bacteria = bact)
  upset(fromList(inpt), order.by = "freq", mainbar.y.label = paste("Core threshold = ",thresh,"%",sep=''))
}
dev.off()