
invsq <- function(i, r1, r2)
{
	(r1[i]-r2[i])/(sum((r1-r2)^2))^(3/2)
}
cross <- function(r1, r2)
{
	c(r1[2]*r2[3]-r1[3]*r2[2], r1[3]*r2[1]-r1[1]*r2[3], r1[1]*r2[2]-r1[2]*r2[1])
}

angles <- function(x, y, z)
{
	u <- diff(x)
	v <- diff(y)
	w <- diff(z)
	u1 <- u[1:(length(u)-1)]
	v1 <- u[1:(length(v)-1)]
	w1 <- u[1:(length(w)-1)]
	u2 <- u[2:length(u)]
	v2 <- v[2:length(v)]
	w2 <- w[2:length(w)]
	theta <- (u1*u2+v1*v2+w1*w2)/sqrt(u1^2+v1^2+w1^2)/sqrt(u2^2+v2^2+w2^2)
}

state <- function(traj, time)
{
	ind <- which.min(abs(traj$time-time))
	c( 
	  traj$X1[ind], traj$X2[ind], traj$X3[ind], 
	  traj$X4[ind], traj$X5[ind], traj$X6[ind]
	  )
}

thrust <- function(parms)
{
	c(0,0,0)
#	with(as.list(parms), {
#		0.01*r0dot/sqrt(sum(r0dot^2))
#	})
}

body3 <- function(t, y, parms, ...)
{
	# cat(y, "\n")
	with( as.list(parms), {
	     	r0 <- c(y[1], y[2], y[3])
		r0dot <- c(y[4], y[5], y[6])

		udot <- y[4]
		vdot <- y[5]
		wdot <- y[6]
		parms <- list(parms, r0=r0, r0dot=r0dot, udot=udot, vdot=vdot, wdot=wdot)

		uddot <- -G*m1*invsq(1, r0, r1)-G*m2*invsq(1, r0, r2) +(2*cross(r0dot, omega) +cross(cross(omega, r0), omega))[1]
		vddot <- -G*m1*invsq(2, r0, r1)-G*m2*invsq(2, r0, r2) +(2*cross(r0dot, omega) +cross(cross(omega, r0), omega))[2]
		wddot <- -G*m1*invsq(3, r0, r1)-G*m2*invsq(3, r0, r2) +(2*cross(r0dot, omega) +cross(cross(omega, r0), omega))[3]
		# cat("F_CF_x=", cross(cross(omega, r0), omega)[1], "\n")
		# cat("F_CF_y=", cross(cross(omega, r0), omega)[2], "\n")
		# cat("F_CF_z=", cross(cross(omega, r0), omega)[3], "\n")
		# cat("F_Cor_x", 2*cross(r0dot, omega)[1],"\n")
		# cat("F_Cor_y", 2*cross(r0dot, omega)[2],"\n")
		# cat("F_Cor_z", 2*cross(r0dot, omega)[3],"\n")
		# cat("F_m2_x", -G*m2*invsq(1, r0, r2),"\n")
		# cat("F_m2_y", -G*m2*invsq(2, r0, r2),"\n")
		# cat("F_m2_z", -G*m2*invsq(3, r0, r2),"\n")
		# cat("F_m1_x", -G*m1*invsq(1, r0, r1),"\n")
		# cat("F_m1_y", -G*m1*invsq(2, r0, r1),"\n")
		# cat("F_m1_z", -G*m1*invsq(3, r0, r1),"\n")
		# cat("omega=", omega, "\n")
		# stop("break")
		# cat("uddot=", uddot, "\nvddot=", vddot, "\nwddot=", wddot, "\n")
		list(c(udot, vdot, wdot, uddot, vddot, wddot))
	})
}

#### Definte the potential ####
pot.grav <- function(x,y,z, x0,y0,z0, M) -M*G/sqrt((x0-x)^2+(y0-y)^2+(z0-z)^2)
pot.CF <- function(x, y, omega.z) -1/2*omega.z^2*((x)^2+(y)^2)
potential <- function(x,y,z) pot.CF(x, y, omega.z)+pot.grav(x,y,z, R1, 0, 0, M1)+pot.grav(x,y,z, R2, 0, 0, M2)
potential.cap <- function(x,y, pot.floor)
{
	cat("Inside potential.cap()\n")
	pot <- potential(x,y,0)
	pot <- pot*(pot >= pot.floor ) +pot.floor*(pot < pot.floor)
}

rainbowPlot <- function()
{
	par.old <- par()
	#### Calculate potential field and make a color key ####
	x <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	y <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	# value below which the potential is not drawn:
	pot.floor <- potential(max(x), max(y), 0)
	cat("\tCalculating potential field.\n")
	pot.XY <- outer(x, y, function(x,y) potential.cap(x, y, pot.floor))
	cols <- seq(from= min(pot.XY, na.rm=T), to= max(pot.XY, na.rm=T), length.out= n.colors)
	cols.mat <- matrix(cols, nrow=1)

	#### Plot the potential field ###
	cat("\tPlotting the potential field.\n")
	nf <- layout(matrix(c(1,2), nrow=1), c(5, 1) )
	image(
	      x, y, pot.XY, 
	      asp=1, col= pal, 
	      bty="n", axes= F, xlab=NA, ylab=NA,
	      # main="Synodic Trajectory"
	      mai="Synodic Coordinates\n(rotational reference frame)"
	      )
	E0 <- 1/2*sum(yini[4:6]^2) +potential(yini[1], yini[2], yini[3])
	if (E0 > potential(x.L4, y.L4, 0))
		cat("\tThere are no bounded regions.\n")
	else
	{
		cat("\tThere are bounded regions.\n")
		contour(x, y, pot.XY, col="black", drawlabels=F, add=T, levels=E0)
	}
	# contour(x, y, pot.XY, add=T, drawlabels=F, levels=c(-1, -1.5, -1.70, -2))
	points(
	       c(0,R1,R2,x.L4,x.L5,x.L1,x.L2,x.L3), 
	       c(0,0,0,y.L4,y.L5,0,0,0), 
	       pch=c(8, 1, 10, 8, 8, 8,8,8)
	       )
	points(R1, 0, pch=".")
	text(
	     x=c(0,R1,R2,x.L4,x.L5,x.L1,x.L2,x.L3), 
	     y=c(0,0,0,y.L4,y.L5,0,0,0), 
	     labels=c("CoM", "M1", "M2", "L4", "L5", "L1", "L2", "L3"), 
	     pos=3, offset=0.3, cex=0.8, font=2
	     )

	#### plot the trajectory of the satellite ####
	cat("\tPlotting the trajectory.\n")
	# plot(traj$X1, traj$X2, asp=1, main= "trjectory plot", xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), type="l", col=alpha("black", 0.5), add=T)
	lines(traj$X1, traj$X2, asp=1, col="lightgray")
	# points(traj$X1, traj$X2, col=alpha("red", 0.01), pch=".")
	points(traj$X1, traj$X2, col="darkred", pch=".")
	# points(rep(parms$r1[1], 2), rep(parms$r1[2], 2), pch=c(1, 19))
	# points(parms$r1[1],parms$r1[2], pch=1)
	# points(parms$r1[1],parms$r1[2], pch=".")
	# points(parms$r2[1],parms$r2[2], pch=10)

	#### plot the color key ####
	cat("\tPlotting the color key.\n")
	par(mar=c(4,0,0,4))
	image(
	      1, cols, cols.mat, 
	      col=pal, bty="n", axes= F,
	      xlab= "U(r)", zlim= c(pot.floor, max(pot.XY, na.rm=T))
	      )
	axis(4, at=c(min(cols), mean(cols), max(cols)), labels= c("LO", "MID", "HI"))

	layout(matrix(1))
	suppressWarnings(par(par.old))
}

# trueTraj <- function(traj, omega.z)
# {
# 	for (i in 1:length(traj$time))
# 	{
# 		t <- traj$time[i]
# 		traj$x.true[i] <- traj$X1[i]*cos(omega.z*t) -traj$X2[i]*sin(omega.z*t)
# 		traj$y.true[i] <- traj$X1[i]*sin(omega.z*t) +traj$X2[i]*cos(omega.z*t)
# 		traj$t[i] <- t
# 	}
# 	return(traj)
# }
trueTraj <- function(traj, omega.z)
{
	x.true <- traj$X1*cos(omega.z*traj$time) -traj$X2*sin(omega.z*traj$time)
	y.true <- traj$X1*sin(omega.z*traj$time) +traj$X2*cos(omega.z*traj$time)
	return(data.frame(x.true, y.true, t=traj$time))
}

Eoft <- function(t, traj)
{
	i <- which.min(traj$time-t)
	1/2*(traj$X4[i]^2+traj$X5[i]^2+traj$X6[i]^2) +potential(traj$X1[i], traj$X2[i], traj$X3[i])
}

library("animation")
animateOrbit <- function(data.in, mp4.out, runtime=10, FPS=20, startTime=0, stopTime= NA)
{
	attach(data.in, name="dat")
	nmax <- runtime*FPS
	interval <- 1/FPS
	if (is.na(stopTime))
		stopTime <- max(traj$time)
	i.start <- which.min(abs(startTime-traj$time))
	i.end <- which.min(abs(stopTime-traj$time))
	ind <- round(seq(from=i.start, to=i.end, length.out=nmax))
	cat("range of ind:\n",range(ind), "\n")
	#### Calculate potential field and make a color key ####
	x <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	y <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	# value below which the potential is not drawn:
	pot.floor <- potential(max(x), max(y), 0)
	cat("\tCalculating potential field.\n")
	pot.XY <- outer(x, y, function(x,y) potential.cap(x, y, pot.floor))
	cols <- seq(from= min(pot.XY, na.rm=T), to= max(pot.XY, na.rm=T), length.out= n.colors)
	cols.mat <- matrix(cols, nrow=1)
	E0 <- 1/2*sum(yini[4:6]^2) +potential(yini[1], yini[2], yini[3])
	if (E0 > potential(x.L4, y.L4, 0))
		cat("\tThere are no bounded regions.\n")
	else
		cat("\tThere are bounded regions.\n")
#	x.L4 <- (R2-R1)*cos(pi/3)+R1
#	y.L4 <- (R2-R1)*sin(pi/3)
#	x.L5 <- x.L4
#	y.L5 <- -y.L4
#	xL1 <- seq(from=R1, to= R2, length.out=1e3)
#	pot1 <- potential(xL1, 0,0)
#	x.L1 <- xL1[which.max(pot1)]
#	xL2 <- seq(from=R2, to= R2+R, length.out=1e3)
#	pot2 <- potential(xL2, 0,0)
#	x.L2 <- xL2[which.max(pot2)]
#	xL3 <- seq(from=R1-R, to= R1, length.out=1e3)
#	pot3 <- potential(xL3, 0,0)
#	x.L3 <- xL3[which.max(pot3)]

	#### Plot the potential field ###
 	saveVideo({
#	saveGIF({
		par.old <- par(mar=c(0,0,0,0))
		ani.options(interval = interval, nmax = nmax)
		for (i in ind)
		{
			image(
			      x, y, pot.XY, 
			      asp=1, col= pal, 
			      bty="n", axes= F, xlab=NA, ylab=NA,
			      # main="Synodic Trajectory"
			      mai="Synodic Coordinates\n(rotational reference frame)"
			      )
			contour(x, y, pot.XY, col="black", drawlabels=F, add=T, levels=E0)
			points(
			       c(0,R1,R2,x.L4,x.L5,x.L1,x.L2,x.L3), 
			       c(0,0,0,y.L4,y.L5,0,0,0), 
			       pch=c(8, 1, 10, 8, 8, 8,8,8)
			       )
			points(R1, 0, pch=".")
			text(
			     x=c(0,R1,R2,x.L4,x.L5,x.L1,x.L2,x.L3), 
			     y=c(0,0,0,y.L4,y.L5,0,0,0), 
			     labels=c("CoM", "M1", "M2", "L4", "L5", "L1", "L2", "L3"), 
			     pos=3, offset=0.3, cex=0.8, font=2
			     )
			tail <- seq(from=abs(i-1e3), to=i)
			lines(traj$X1[tail], traj$X2[tail], asp=1, col="black", lwd=2)
			lines(traj$X1[tail], traj$X2[tail], asp=1, col="white", lwd=0.7)
			points(traj$X1[i], traj$X2[i], col="darkred", bg="white", pch=21, cex=3)
		}
		par(par.old)
	}, video.name=mp4.out, other.opts = "-pix_fmt yuv420p -b 300k")
#	}, movie.name= mp4.out)
	detach(dat, unload=T)
}

