strain =c(
374/4699,
438/5055,
328/3220,
353/4216,
463/4518,
468/4515,
512/4374,
389/4468,
337/4199,
380/4094,
447/4317,
520/4399,
447/4494,
298/4149
)


speciesInt = c( 
325/4564,
328/3220,
353/4216,
289/3588,
253/3461
)


speciesAII = c( 
328/3220,
353/4216,
346/4564,
319/3588,
337/3601
)

speciesDollo = c( 
328/3220,
353/4216,
325/4564,
354/3892,
365/3959
)


subFamilyInt = c(
280/2960,
233/2604,
251/2668
)

subFamilyAII = c(
323/2960,
327/2605,
274/2668

)

subFamilyDollo = c(
336/3408,
409/4120,
284/2856
)

familyInt = c(
184/1908
)

familyAII = c(
273/1908
)

familyDollo = c(
410/4344
)


pdf(file = '~/EnTrI/figures/method-comparison.pdf')

plot(4*strain/strain, strain, xlim=c(1,4), ylim=c(0.07,0.15), col="olivedrab4", ylab="proportion of ancestral genes", xaxt="n", xlab="",pch=19)

points(3*speciesInt/speciesInt,speciesInt, col="red4", pch=19)
points(3*speciesAII/speciesAII,speciesAII, col="olivedrab4", pch=19)
points(3*speciesDollo/speciesDollo,speciesDollo, col="cyan3", pch=19)


points(2*subFamilyInt/subFamilyInt,subFamilyInt, col="red4", pch=19)
points(2*subFamilyAII/subFamilyAII,subFamilyAII, col="olivedrab4", pch=19)
points(2*subFamilyDollo/subFamilyDollo,subFamilyDollo, col="cyan3", pch=19)

points(1*familyInt/familyInt,familyInt, col="red4", pch=19)
points(1*familyAII/familyAII,familyAII, col="olivedrab4", pch=19)
points(1*familyDollo/familyDollo,familyDollo, col="cyan3", pch=19)

lines(1:4, c(median(familyInt), median(subFamilyInt), median(speciesInt), median(strain) ), col="red4", lwd=3)

lines(1:4, c(median(familyAII), median(subFamilyAII), median(speciesAII), median(strain) ), col="olivedrab4", lwd=3)
lines(1:4, c(median(familyDollo), median(subFamilyDollo), median(speciesDollo), median(strain) ), col="cyan3", lwd=3)

legend("topright", c("intersection", "ancestral II", "Dollo"), fil=c("red4", "olivedrab4", "cyan3"))
axis(1, at = 1:4, labels=c("family", "subfam.", "species", "strain"))

dev.off()