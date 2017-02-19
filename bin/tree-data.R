strain =c(
  357/4571,
  315/4988,
  373/5672,
  325/4957,
  314/4520,
  303/4517,
  316/4492,
  307/4453,
  302/4312,
  308/4322,
  308/4980,
  321/4805,
  327/4621,
  345/4735,
  296/4306,
  299/4317
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
357/4571,
287/4563,
325/4957,
269/3843,
273/3829
)


subFamilyFitch = c(
268/3251,
280/3339,
279/3411
)


familyFitch = c(
288/3169
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
