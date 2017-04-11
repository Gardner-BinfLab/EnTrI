datapath <- '~/EnTrI/results/interesting_genes/sometimes-essential-marked-dup-with-modules.tsv'
corr.list <- list()
data <- read.csv(datapath, sep='\t', stringsAsFactors = FALSE)
data <- data[ order(data$path_name, data$Gene), ]
data[data == 'X'] <- -7
# i = 1
# j = i
counter = 0
# while (j <= nrow(data))
# {
#   datachunk <- list()
#   line1 <- data[i,]
#   while (j <= nrow(data) & line1[36] == data[j, 36])
#   {
#     line2 <- data[j,]
#     datachunk <- rbind(datachunk, as.numeric(line2[3:17]))
#     j = j + 1
#   }
#   rownames(datachunk) <- data[i:(j-1),1]
#   # datachunk[datachunk <= 1.644854] <- 0
#   # datachunk[datachunk > 1.644854] <- 1
#   i = j
#   if (nrow(datachunk) > 1)
#   {
#     for (k in seq(1, (nrow(datachunk)-1), 1))
#     {
#       for (l in seq((k+1), nrow(datachunk), 1))
#       {
#         counter = counter + 1
#         correlation <- cor.test(as.numeric(datachunk[k,]),as.numeric(datachunk[l,]), alternative="less", method = "pearson")
#         if (correlation$p.value < 0.05)
#         {
#           corr.list <- rbind(corr.list, c(row.names(datachunk)[k], row.names(datachunk)[l], correlation$p.value))
#         }
#       }
#     }
#   }
# }
for (i in seq(1, nrow(data)-1, 1))
{
     for (j in seq(i, nrow(data), 1))
     {
       counter = counter + 1
       correlation <- cor.test(as.numeric(data[i,3:17]),as.numeric(data[j,3:17]), alternative="less", method = "pearson")
       if (correlation$p.value < 0.05)
       {
         corr.list <- rbind(corr.list, c(data[i,1], data[j,1], correlation$p.value))
       }
     }
}
p.adj <- p.adjust(as.numeric(corr.list[,3]), method = "BH", counter)
corr.list <- cbind(corr.list, p.adj)
corr.list <- as.data.frame(matrix(unlist(corr.list), nrow=nrow(corr.list)))
corr.list <- corr.list[ order(corr.list[,4]), ]
# corr.list.filtered <- unique(corr.list[corr.list[,4]<0.05,])