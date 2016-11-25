
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
pot.grav <- function(x,y, x0, y0, M) -M*G/sqrt((x0-x)^2+(y0-y)^2)
pot.CF <- function(x, y, omega.z) -1/2*omega.z^2*((x)^2+(y)^2)
potential <- function(x,y) pot.CF(x, y, omega.z)+pot.grav(x, y, R1, 0, M1)+pot.grav(x, y, R2, 0, M2)
potential.cap <- function(x,y, pot.floor)
{
	cat("Inside potential.cap()\n")
	pot <- potential(x,y)
	pot <- pot*(pot >= pot.floor ) +pot.floor*(pot < pot.floor)
}

rainbowPlot <- function()
{
	par.old <- par()
	#### Calculate potential field and make a color key ####
	x <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	y <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=3e2)
	# value below which the potential is not drawn:
	pot.floor <- potential(max(x), max(y))
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
	      bty="n", axes= F, xlab=NA, ylab=NA
	      )
	# contour(x, y, pot.XY, col="lightgray", drawlabels=F, add=T)
	points(c(0,R1,R2,x.L4,x.L5), c(0,0,0,y.L4,y.L5), pch=c(8, 1, 10, 8, 8))
	points(R1, 0, pch=".")
	text(x=c(0,R1,R2,x.L4,x.L5), y=c(0,0,0,y.L4,y.L5), labels=c("CoM", "M1", "M2", "L4", "L5"), pos=3, offset=0.3, cex=0.8, font=2)

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
