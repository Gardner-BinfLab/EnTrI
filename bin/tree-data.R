strain =c(
  360/4317,
  303/5003,
  372/5677,
  307/4957,
  335/4525,
  323/4522,
  321/4492,
  371/4455,
  298/4316,
  352/4341,
  311/4981,
  324/4831,
  439/4622,
  377/4764,
  438/4311,
  299/4317
)


speciesFitch = c( 
  360/4317,
  278/4580,
  307/4957,
  302/3852,
  291/3845
)


subFamilyFitch = c(
  278/3159,
  283/3348,
  283/3412
)


familyFitch = c(
  289/3160
)



pdf(file = '~/EnTrI/figures/fitch.pdf')

plot(4*strain/strain, strain, xlim=c(1,4), ylim=c(0.06,0.102), col="red4", ylab="proportion of ancestral genes", xaxt="n", xlab="",pch=19,
     cex.lab=1.5, cex.axis=1.5)

points(3*speciesFitch/speciesFitch,speciesFitch, col="red4", pch=19)

points(2*subFamilyFitch/subFamilyFitch,subFamilyFitch, col="red4", pch=19)

points(1*familyFitch/familyFitch,familyFitch, col="red4", pch=19)

lines(1:4, c(median(familyFitch), median(subFamilyFitch), median(speciesFitch), median(strain) ), col="red4", lwd=3)

#legend("topright", c("intersection", "ancestral II", "Dollo"), fil=c("red4", "olivedrab4", "cyan3"))
axis(1, at = 1:4, labels=c("family", "subfam.", "species", "strain"), cex.axis=1.5)

dev.off()
