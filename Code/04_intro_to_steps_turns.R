## Introduction to Animal Movement Ecology
## K. Manlove, August 2023
##
## Probability models for step length and turning angle

library(geostats)

angles <- seq(0, 2*pi, by = .01)
angles <- seq(from=-pi,to=pi, length.out=200)
x_kappa.1 <- vonMises(angles, mu = pi/4, kappa = .1)
x_kappa01 <- vonMises(angles, mu = pi/4, kappa = 1)
x_kappa05 <- vonMises(angles, mu = pi/4, kappa = 5)
x_kappa10 <- vonMises(angles, mu = pi/4, kappa = 10)

pdf("ExampleScripts/Output/path_segmentation/vonMisesPDF.pdf", height = 4, width = 7)
par(mfrow = c(1, 2))
plot(x_kappa10 ~ angles, type = "l", col = "orange",
     xlab = "radian", ylab = "Pr[X = x]",
     las = 1, lwd = 2)
lines(x_kappa05 ~ angles, type = "l", col = "red", lwd = 2)
lines(x_kappa01 ~ angles, type = "l", col = "purple", lwd = 2)
lines(x_kappa.1 ~ angles, type = "l", col = "blue", lwd = 2)
leg_text <- c("kappa = 0.1", "kappa = 1", "kappa = 5", "kappa = 10")
leg_col <- c("blue", "purple", "red", "orange")
legend("topleft", leg_text, col = leg_col, lwd = rep(2, 4), bty = "n")
mtext(side = 3, outer = F, "Conventional PDF diagram", line = 2)

plot(x=c(-1,1.2), y=c(-1,1.2), type='n',
     axes=FALSE, ann=FALSE, bty='n', asp=1)

#d_kappa01 <- vonMises(a=a, mu=pi/4, kappa=1)
symbols(x=0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
lines(x=(1+x_kappa.1)*cos(a),y=(1+x_kappa.1)*sin(a), xpd=NA, 
      col = "blue", lwd = 2)
lines(x=(1+x_kappa01)*cos(a),y=(1+x_kappa01)*sin(a), xpd=NA, 
      col = "purple", lwd = 2)
lines(x=(1+x_kappa05)*cos(a),y=(1+x_kappa05)*sin(a), xpd=NA, 
      col = "red", lwd = 2)
lines(x=(1+x_kappa10)*cos(a),y=(1+x_kappa10)*sin(a), xpd=NA, 
      col = "orange", lwd = 2)
mtext(side = 3, outer = F, "PDF in polar coordinates", line = 2)
dev.off()
