#!/usr/bin/Rscript
# This is a simulation of the restricted 3-Body problem using the deSolve package.
# See https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Soetaert~et~al.pdf
# for a quick overview of how to use the package, and the official docs for more
# details.
# Joe Shields
# 2016-10-6

# load packages
require("deSolve")
# require("scatterplot3d")
source("orbitFunctions.R")

#### Input parameters ####
# time steps for the simulation:
times <- seq(from=0, by=1e-1, to=1e4) 
# initial conditions c(x0, y0, z0, v0, u0, w0):
yini <- c(
	  25,		44,		0, 
	  0.2,		0,		0.1
	  )
# L1 (needs work)
# yini <- c(
# 	  1.15,		0.1,		0, 
# 	  -1,		1.1,		0
# 	  )
# L3 (needs work)
# yini <- c(
# 	  -2.8,		0,		0, 
# 	  0.1,		0.01,		0
# 	  )
# L4 (needs work)
# yini <- c(
# 	  x.L4,		y.L4+0.1,	0, 
# 	  0.01,		0,		0.01
# 	  )
# L5 (needs work)
# yini <- c(
# 	  x.L5+0.01*runif(1),		y.L5+0.01*runif(1),		0+0.01*runif(1), 
# 	  0+0.01*runif(1),		0+0.01*runif(1),		0+0.01*runif(1)
# 	  )
# strength of gravity:
G <- 10
# mass of the "sun":
M1 <- 100 # 4
# mass of the "planet":
M2 <- 1
# reduced mass for M1 and M2:
mu <- M1*M2/(M1+M2)
# distance from the CoM to M1:
R1 <- -0.5
# distance from the CoM to M2:
R2 <- -M1*R1/M2
# distance between the big masses:
R <- abs(R2-R1)
# rotation rate of the reference frame:
omega.z <- sqrt(G*(M1+M2)/abs(R2-R1)^3)
# how many colors to use to plot the potential:
n.colors <- 2^6
# how much of the violet end of the colors to cut off:
end.colors <- 0.8
# how much of the red end of the colors to cut off:
start.colors <- 0.05
# color palette for the potential:
pal <- rev(rainbow(n.colors, start= start.colors,end=end.colors)) 
# how zoomed-in the plot is (smaller numbers -> smaller view):
zoom.factor <- 2

#### Definte the potential ####
pot.grav <- function(x,y, x0, y0, M) -M*G/sqrt((x0-x)^2+(y0-y)^2)
pot.CF <- function(x, y, omega.z) -1/2*omega.z^2*((x)^2+(y)^2)
potential <- function(x,y) pot.CF(x, y, omega.z)+pot.grav(x, y, R1, 0, M1)+pot.grav(x, y, R2, 0, M2)
potential.cap <- function(x,y)
{
	pot <- potential(x,y)
	pot <- pot*(pot >= pot.floor ) +pot.floor*(pot < pot.floor)
}

#### Calculate potential field and make a color key ####
x <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=5e2)
y <- seq(from=-zoom.factor*R, to=zoom.factor*R, length.out=5e2)
# value below which the potential is not drawn:
pot.floor <- potential(max(x), max(y))
cat("\tCalculating potential field.\n")
pot.XY <- outer(x, y, potential.cap)
cols <- seq(from= min(pot.XY, na.rm=T), to= max(pot.XY, na.rm=T), length.out= n.colors)
cols.mat <- matrix(cols, nrow=1)

#### Calculate the locations of the Lagrange points ####
x.L4 <- (R2-R1)*cos(pi/3)+R1
y.L4 <- (R2-R1)*sin(pi/3)
x.L5 <- x.L4
y.L5 <- -y.L4
# use polyroot() to easily find the locations of L1 and L2

#### Calculate the motion of the satellite ####
parms <- list(r1=c(R1, 0, 0), r2=c(R2, 0, 0), m1=M1, m2=M2, omega= c(0,0,omega.z))
traj <- ode(y= yini, func=body3, times= times, parms=parms, method="rk4")
traj <- data.frame(traj)

#### Check the satellite motion for kinks ####
ang <- asin(angles(traj$X1, traj$X2, traj$X3))*360/2/pi
cat(
    "Maximum step angle", max(ang), "degrees ocurring at step", which.max(ang), 
    "(", traj$time[which.max(ang)],"seconds).\n"
    )

source("orbitPlot.R")
# system("cp orbit.pdf ~/Downloads/")
