strain =c(
  307/4957,
  352/4341,
  298/4316,
  335/4525,
  323/4522,
  321/4492,
  371/4455,
  340/4831,
  334/4981,
  438/4311,
  299/4317,
  388/5677,
  303/5003,
  360/4317
)


speciesFitch = c( 
  307/4957,
  302/3853,
  309/3708,
  277/4524,
  360/4317
)


subFamilyFitch = c(
  283/3349,
  309/3708,
  277/3084
)


familyFitch = c(
  296/3164
)



pdf(file = '~/EnTrI/figures/fitch.pdf')

par(mar=c(5.1,5.1,4.1,2.1))

plot(4*strain/strain, strain, xlim=c(1,4), ylim=c(0.06,0.102), col="red4", ylab="Proportion of essential ancestral genes", xaxt="n", xlab="",pch=19,
     cex.lab=2, cex.axis=2)

points(3*speciesFitch/speciesFitch,speciesFitch, col="red4", pch=19)

points(2*subFamilyFitch/subFamilyFitch,subFamilyFitch, col="red4", pch=19)

points(1*familyFitch/familyFitch,familyFitch, col="red4", pch=19)

lines(1:4, c(median(familyFitch), median(subFamilyFitch), median(speciesFitch), median(strain) ), col="red4", lwd=3)

#legend("topright", c("intersection", "ancestral II", "Dollo"), fil=c("red4", "olivedrab4", "cyan3"))
axis(1, at = 1:4, labels=c("family", "subfam.", "species", "strain"), cex.axis=2)

dev.off()
