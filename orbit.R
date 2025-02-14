#!/usr/bin/Rscript
# This is a simulation of the restricted 3-Body problem using the deSolve package.
# See https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Soetaert~et~al.pdf
# for a quick overview of how to use the package, and the official docs for more
# details.
# Use something similar to the following to render videos.
# animateOrbit(data.in="interestingOrbits/archive/Tue_Nov_29_23-21-27_2016.RData", mp4.out="confinedTransfers4.gif", runtime=60, FPS=30)
# Joe Shields
# 2016-10-6

# load packages
require("deSolve")
# require("scatterplot3d")
source("orbitFunctions.R")

#### Input parameters ####
# time steps for the simulation:
# Use system.time(source("orbit.R")) when moving to large simulation times 
# to estimate run time as a function of steps.
times <- seq(from=0, by=0.0001, to=60) 
# initial conditions c(x0, y0, z0, v0, u0, w0):
yini <- c(
	  0.2865,	0.0,		0, 
	  0.41,		-0.39,		0
	  )
mu.2 <- 0.3
# Normalized parameters:
M2 <- mu.2
M1 <- 1-mu.2
R1 <- -M2
R2 <- M1
R <- 1
# # strength of gravity:
G <- 1
# # mass of the "sun":
# M1 <- 100 # 4
# # mass of the "planet":
# M2 <- 1
# # reduced mass for M1 and M2:
# mu <- M1*M2/(M1+M2)
# # distance from the CoM to M1:
# R1 <- -0.5
# # distance from the CoM to M2:
# R2 <- -M1*R1/M2
# # distance between the big masses:
# R <- abs(R2-R1)
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
zoom.factor <- 1.5

#### Calculate the locations of the Lagrange points ####
x.L4 <- (R2-R1)*cos(pi/3)+R1
y.L4 <- (R2-R1)*sin(pi/3)
x.L5 <- x.L4
y.L5 <- -y.L4
xL1 <- seq(from=R1, to= R2, length.out=1e3)
pot1 <- potential(xL1, 0,0)
x.L1 <- xL1[which.max(pot1)]
xL2 <- seq(from=R2, to= R2+R, length.out=1e3)
pot2 <- potential(xL2, 0,0)
x.L2 <- xL2[which.max(pot2)]
xL3 <- seq(from=R1-R, to= R1, length.out=1e3)
pot3 <- potential(xL3, 0,0)
x.L3 <- xL3[which.max(pot3)]
# use polyroot() to easily find the locations of L1 and L2

#### Calculate the motion of the satellite ####
cat("\tCalculating the motion of the satellite. This may take a while. (", length(times),"steps)\n")
cat("\tEstimated crunch time (based on large runs on my laptop):", length(times)*712e-6, "seconds\n")
parms <- list(r1=c(R1, 0, 0), r2=c(R2, 0, 0), m1=M1, m2=M2, omega= c(0,0,omega.z))
traj <- ode(y= yini, func=body3, times= times, parms=parms, method="rk4")
traj <- data.frame(traj)

#### Check the satellite motion for kinks ####
ang <- asin(angles(traj$X1, traj$X2, traj$X3))*360/2/pi
cat(
    "Maximum step angle", max(ang), "degrees ocurring at step", which.max(ang), 
    "(", traj$time[which.max(ang)],"seconds).\n"
    )

dateStr <- gsub(" ", "_", date())
dateStr <- paste("interestingOrbits/archive/", dateStr, ".RData", sep="")
dateStr <- gsub(":", "-", dateStr)
cat("\tSaving simulation to", dateStr, "\n")
save.image(file= dateStr)
# Use load() to recover and plot these. 

source("orbitPlot.R")
# system("cp orbit.pdf ~/Downloads/")
