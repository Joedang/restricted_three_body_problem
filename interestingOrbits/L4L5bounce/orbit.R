#!/usr/bin/Rscript
if (require("deSolve")!=T)
	install.packages("deSolve")
pdf("orbit.pdf")
# DHM <- function(t, y, parms, ...)
# {
# 	with( as.list(parms), {
# 		list(c(
# 		       y[2], 
# 		       -k1*y[1]-mu1*y[2], 
# 		       y[4], 
# 		       -k2*y[3]-mu1*y[4]
# 		       ))
# 	})
# }
# yini <- c(1, 0, 1, 1)
# parameters <- c(mu1=0.1, mu2= 1, k1=2, k2=2)
# vib <- ode(y=yini, func=DHM, parms= parameters, times= seq(from=0, to=10, length.out=1e3))
# vib <- data.frame(vib)
# # plot(vib$time, vib$X1)
# # plot(vib$time, vib$X3)
# plot(vib$X2, vib$X3, main="config space")
# 
# # This actually works! (yay)
# # dev.off()
# 
G <- 1
M <- 1e1
# orbit <- function(t, y, ...)
# {
# 	udot <- y[2]
# 	uddot <- -G*M*y[1]/(y[1]^2+y[3]^2)^(3/2)
# 	vdot <- y[4]
# 	vddot <- -G*M*y[3]/(y[1]^2+y[3]^2)^(3/2)
# 	list(c(udot, uddot, vdot, vddot))
# }
# yini <- c(1,0,0,2)
# path <- ode(y= yini, func=orbit, times= seq(from=0, to=5, length.out=1e3), parms=1, method="rk4")
# path <- data.frame(path)
# # plot(path$time, path$X1)
# # plot(path$time, path$X3)
# plot(path$X1, path$X3, asp=1)
# points(0,0,pch=23)

invsq <- function(i, r1, r2)
{
	(r1[i]-r2[i])/(sum((r1-r2)^2))^(3/2)
}
cross <- function(r1, r2)
{
	c(r1[2]*r2[3]-r1[3]*r2[2], r1[3]*r2[1]-r1[1]*r2[3], r1[1]*r2[2]-r1[2]*r2[1])
}

body3 <- function(t, y, parms, ...)
{
	with( as.list(parms), {
	     	r0 <- c(y[1], y[2], 0)
		r0dot <- c(y[3], y[4], 0)

		udot <- y[3]
		vdot <- y[4]
		uddot <- -G*m1*invsq(1, r0, r1)-G*m2*invsq(1, r0, r2) +(2*cross(r0dot, omega) +cross(cross(omega, r0), omega))[1]
		vddot <- -G*m1*invsq(2, r0, r1)-G*m2*invsq(2, r0, r2) +(2*cross(r0dot, omega) +cross(cross(omega, r0), omega))[2]
		list(c(udot, vdot, uddot, vddot))
	})
}
yini <- c(cos(pi/3), sin(pi/3), 0.01, 0)
parms <- list(r1=c(0, 0, 0), r2=c(1, 0, 0), m1=10, m2=0.1, omega= c(0,0,1)*sqrt(G*10/1))
traj <- ode(y= yini, func=body3, times= seq(from=0, by=5e-3, length.out=10e3), parms=parms, method="rk4")
traj <- data.frame(traj)
# plot(traj$time, traj$X1)
# plot(traj$time, traj$X2)
plot(traj$X1, traj$X2, asp=1, main= "trjectory plot", xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), type="l")
points(parms$r1[1],parms$r1[2], pch=23)
points(parms$r2[1],parms$r2[2], pch=23)
# lines(path$X1, path$X3)

dev.off()
system("cp orbit.pdf ~/Downloads/")
