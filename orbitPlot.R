library("animation")
if (!exists("GIF"))
	GIF <- F
# orbitPlot <- function()
# {
pdf("orbit.pdf")
par.def <- par(no.readonly=T)

plot(traj$X1, traj$X2, asp=1, main= "trjectory plot", xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), type="l")
# points(rep(parms$r1[1], 2), rep(parms$r1[2], 2), pch=c(1, 19))
points(parms$r1[1],parms$r1[2], pch=1)
points(parms$r1[1],parms$r1[2], pch=".")
points(parms$r2[1],parms$r2[2], pch=10)

# time domain of XYZ coords in rotating frame
layout(matrix(1:3, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(traj$time, traj$X1, type="l", col="red", ylab="x")
par(mar=c(0,4,0,2))
plot(traj$time, traj$X2, type="l", col="green", ylab="y")
par(mar=c(5,4,0,2))
plot(traj$time, traj$X3, type="l", col="blue", ylab="z", xlab="t")
par(par.old)

layout(matrix(1:3, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(density(traj$X1))
plot(density(traj$X2))
plot(density(traj$X3))
par(par.old)

layout(matrix(1:2, ncol=1))
par.old <- par(mar=c(0,4,4,2))
v <- c(0, sqrt(diff(traj$X1)^2 +diff(traj$X2)^2 +diff(traj$X3)^2 ))/diff(traj$time)[1]
KE <- 0.5*v^2
plot(traj$time, KE, main="Kinetic Energy", type="l", ylab="KE/m")
plot(density(KE))
par(par.old)

layout(matrix(1:2, ncol=1))
par.old <- par(mar=c(0,4,4,2))
plot(traj$X1, traj$X3, asp=1, type="l")
plot(traj$X2, traj$X3, asp=1, type="l")
par(par.old)

par(par.def)
scatterplot3d(traj$X1, traj$X2, traj$X3, scale.y=2, type="l")

dev.off()

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
# system("")

system("cp orbit.pdf ~/Downloads/")
# }
