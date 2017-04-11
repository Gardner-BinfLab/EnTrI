strain =c(
  296/3562,
  316/4719,
  359/5081,
  317/4236,
  314/4523,
  303/4520,
  314/4468,
  303/4374,
  298/4202,
  303/4109,
  309/4783,
  327/4739,
  318/4407,
  331/4530,
  291/4290,
  299/4300
)


speciesFitch = c( 
  296/3562,
  293/4580,
  317/4236,
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
