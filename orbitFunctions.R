
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
