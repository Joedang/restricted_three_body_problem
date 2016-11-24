# library("animation")
require("scales")
if (!exists("GIF"))
{
	GIF <- F
       	cat("Defaulting to no GIF creation...\n")
}
pdf("orbit.pdf")
par.def <- par(no.readonly=T)

#### Plot the potential field ###
nf <- layout(matrix(c(1,2), nrow=1), c(5, 1) )
image(
      x, y, pot.XY, 
      asp=1, col= pal, 
      bty="n", axes= F, xlab=NA, ylab=NA
      )
# contour(x, y, pot.XY, col="lightgray", drawlabels=F, add=T)
points(c(0,R1,R2,x.L4,x.L5), c(0,0,0,y.L4,y.L5), pch=c(8, 8, 8))
text(x=c(0,R1,R2,x.L4,x.L5), y=c(0,0,0,y.L4,y.L5), labels=c("CoM", "M1", "M2", "L4", "L5"), pos=3, offset=0.3, cex=0.8, font=2)

#### plot the trajectory of the satellite ####
plot(traj$X1, traj$X2, asp=1, main= "trjectory plot", xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), type="l", col=alpha("black", 0.5), add=T)
points(traj$X1, traj$X2, col=alpha("red", 0.01), pch=".")
# points(rep(parms$r1[1], 2), rep(parms$r1[2], 2), pch=c(1, 19))
points(parms$r1[1],parms$r1[2], pch=1)
points(parms$r1[1],parms$r1[2], pch=".")
points(parms$r2[1],parms$r2[2], pch=10)

#### plot the color key ####
par.old <- par(mar=c(5,0,5,3))
image(
      1, cols, cols.mat, 
      col=pal, bty="n", axes= F,
      xlab= "U(r)", zlim= c(pot.floor, max(pot.XY, na.rm=T))
      )
axis(4, at=c(min(cols), mean(cols), max(cols)), labels= c("LO", "MID", "HI"))
par(par.old)


#### time domain of XYZ coords in rotating frame ####
layout(matrix(1:3, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(traj$time, traj$X1, type="l", col="red", ylab="x")
par(mar=c(0,4,0,2))
plot(traj$time, traj$X2, type="l", col="green", ylab="y")
par(mar=c(5,4,0,2))
plot(traj$time, traj$X3, type="l", col="blue", ylab="z", xlab="t")
par(par.old)

#### smoothed histogram of X, Y, and Z coordinates ####
layout(matrix(1:3, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(density(traj$X1))
plot(density(traj$X2))
plot(density(traj$X3))
par(par.old)

#### time domain and smoothed histogram of the kinetic energy ####
layout(matrix(1:2, ncol=1))
par.old <- par(mar=c(0,4,4,2))
v <- c(0, sqrt(diff(traj$X1)^2 +diff(traj$X2)^2 +diff(traj$X3)^2 ))/diff(traj$time)[1]
KE <- 0.5*v^2
plot(traj$time, KE, main="Kinetic Energy", type="l", ylab="KE/m")
plot(density(KE))
par(par.old)

#### Z as a function of X, also as a function of Y ####
layout(matrix(1:2, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(traj$X1, traj$X3, asp=1, type="l")
plot(traj$X2, traj$X3, asp=1, type="l")
par(par.old)

par(par.def)
scatterplot3d(traj$X1, traj$X2, traj$X3, scale.y=2, type="l")

dev.off()

#### make a GIF of the motion ####
if (GIF)
{
	pdf("/tmp/orb%03d.pdf", onefile=F)
	for (i in round(seq(from=1, to=length(traj$X1), length.out=200)))
	{
		s3d <- scatterplot3d(
				     traj$X1[1:i],
				     traj$X2[1:i], 
				     traj$X3[1:i],
				     scale.y=2, type="l", angle=80, 
				     xlim=range(traj$X1), 
				     ylim=range(traj$X2), 
				     zlim=range(traj$X3), 
				     )
		s3d$points3d(traj$X1[i], traj$X2[i], traj$X3[i])
	}
	dev.off()
	system("bash orbitGIF.sh")
}

# system("cp orbit.pdf ~/Downloads/")
