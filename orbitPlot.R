# library("animation")
# require("scales")
cat("\tBegining the plotting sequence.\n")
if (!exists("GIF"))
{
	GIF <- F
       	cat("\tDefaulting to no GIF creation...\n")
}
pdf("orbit.pdf")
par.def <- par(no.readonly=T)

#### lay down a plot of the trajectory with labeled points and potential field ####
rainbowPlot()

#### time domain of XYZ coords in rotating frame ####
cat("\tPlotting the other heuristics.\n")
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

ttj <- trueTraj(traj, omega.z)
plot(ttj$x.true, ttj$y.true, asp=1, type="l")
# scatterplot3d(traj$X1, traj$X2, traj$X3, scale.y=2, type="l")

#### report initial values ####
cat("\tReporting initial conditions.\n")
plot.new()
# report <- NULL
# for (vari in ls())
# {
# 	val=get(vari)
# 	if (typeof(val)!="closure")
# 	{
# 		if (length(val) > 6 )
# 			val=head(val, 6)
# 		# report <- paste(report, vari, "=", val, "\n")
# 		cat(paste(report, vari, "=", val, "\n"))
# 	}
# }
# cat(report)
text(0.2,0.5,paste(
		   "yini=", capture.output(print(yini)), 
		   "\nM1=", M1, 
		   "\nM2=", M2, 
		   "\nG=", G, 
		   "\nR1=", R1, 
		   "\ndateStr=", dateStr
		   ), pos=4)

dev.off()

#### make a GIF of the motion ####
if (GIF)
{
	cat("\tWriting a GIF. This usually takes a while.\n")
	pdf("/tmp/orb%03d.pdf", onefile=F)
	for (i in round(seq(from=1, to=length(traj$X1), length.out=200)))
	{
		# not ready yet; re-calculates potential field every time; makes PDFs take a long time to convert to JPG
		# I really need to just recompile R to allow Yihui's GIF library.
		# rainbowPlot() 
		# s3d <- scatterplot3d(
		# 		     traj$X1[1:i],
		# 		     traj$X2[1:i], 
		# 		     traj$X3[1:i],
		# 		     scale.y=2, type="l", angle=80, 
		# 		     xlim=range(traj$X1), 
		# 		     ylim=range(traj$X2), 
		# 		     zlim=range(traj$X3), 
		# 		     )
		# s3d$points3d(traj$X1[i], traj$X2[i], traj$X3[i])
	}
	dev.off()
	system("bash orbitGIF.sh")
}

cat("\tCopying output to ~/Downloads.\n")
system("cp orbit.pdf ~/Downloads/")
