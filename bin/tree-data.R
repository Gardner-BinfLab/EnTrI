strain =c(
  306/3713,
  315/4761,
  359/5113,
  316/4239,
  314/4518,
  304/4515,
  315/4468,
  302/4374,
  301/4195,
  301/4096,
  309/4783,
  327/4739,
  322/4461,
  333/4508,
  289/4289,
  299/4300
)

# strain =c(
# 318/4988,
# 290/4306,
# 327/4621,
# 311/4980,
# 282/3725,
# 370/5672,
# 344/4735,
# 325/4805,
# 329/4957,
# 303/4312,
# 314/4520,
# 304/4517,
# 317/4492,
# 301/4453,
# 306/4322
# )


speciesFitch = c( 
  306/3713,
  291/4618,
  316/4239,
  268/3843,
  279/3870
)


subFamilyFitch = c(
  280/3339,
  272/3271,
  285/3339
)


familyFitch = c(
  291/3178
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
