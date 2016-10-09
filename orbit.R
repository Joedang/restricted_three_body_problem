#!/usr/bin/Rscript
# This is a simulation of the restricted 3-Body problem using the deSolve package.
# See https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Soetaert~et~al.pdf
# for a quick overview of how to use the package, and the official docs for more
# details.
# Joe Shields
# 2016-10-6

# load packages
require("deSolve")
require("scatterplot3d")

G <- 1
source("orbitFunctions.R")
#### 3-Body Problem ####
# bounce between L4 and L5:
# yini <- c(cos(pi/3), sin(pi/3), 0, 0.01*runif(1), 0.01*runif(1), 0.01*runif(1))
# parms <- list(r1=c(0, 0, 0), r2=c(1, 0, 0), m1=10, m2=0.1, omega= c(0,0,1)*sqrt(G*10/1))
# tiny L3 orbit:
# yini <- c(
# 	  -1+7*0.1/12/10,	0.1,		0, 
# 	  0,	0,	0.1
# 	  )
# parms <- list(r1=c(0, 0, 0), r2=c(1, 0, 0), m1=10, m2=0.1, omega= c(0,0,1)*sqrt(G*10/1))
# yini <- c(0.6, 0, 0, 0.01*runif(1), 2+0.01*runif(1), 0.01*runif(1))
parms <- list(r1=c(0, 0, 0), r2=c(1, 0, 0), m1=10, m2=0.1, omega= c(0,0,1)*sqrt(G*10/1))
traj <- ode(y= yini, func=body3, times= seq(from=0, by=2e-3, to=3.87), parms=parms, method="rk4")
traj <- data.frame(traj)

ang <- asin(angles(traj$X1, traj$X2, traj$X3))*360/2/pi
cat(
    "Maximum step angle", max(ang), "degrees ocurring at step", which.max(ang), 
    "(", traj$time[which.max(ang)],"seconds).\n"
    )

source("orbitPlot.R")
# orbitPlot()
