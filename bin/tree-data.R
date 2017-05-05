strain =c(
  336/4317,
  319/5003,
  372/5677,
  330/4957,
  314/4525,
  303/4522,
  315/4492,
  307/4455,
  300/4316,
  310/4341,
  311/4981,
  324/4831,
  328/4622,
  344/4764,
  291/4311,
  299/4317
)


speciesFitch = c( 
  336/4317,
  293/4580,
  330/4957,
  267/3852,
  275/3845
)


subFamilyFitch = c(
  265/3159,
  285/3348,
  283/3412
)


familyFitch = c(
  287/3160
)



pdf(file = '~/EnTrI/figures/fitch.pdf')

plot(4*strain/strain, strain, xlim=c(1,4), ylim=c(0.06,0.092), col="red4", ylab="proportion of ancestral genes", xaxt="n", xlab="",pch=19,
     cex.lab=1.5, cex.axis=1.5)

points(3*speciesFitch/speciesFitch,speciesFitch, col="red4", pch=19)

points(2*subFamilyFitch/subFamilyFitch,subFamilyFitch, col="red4", pch=19)

points(1*familyFitch/familyFitch,familyFitch, col="red4", pch=19)

lines(1:4, c(median(familyFitch), median(subFamilyFitch), median(speciesFitch), median(strain) ), col="red4", lwd=3)

#legend("topright", c("intersection", "ancestral II", "Dollo"), fil=c("red4", "olivedrab4", "cyan3"))
axis(1, at = 1:4, labels=c("family", "subfam.", "species", "strain"), cex.axis=1.5)

dev.off()
