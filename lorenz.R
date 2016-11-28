require("deSolve")

lorenz <- function(t, y, parms)
{
	with( as.list(parms), {
	     xdot <- sigma*(y[2]-y[1])
	     ydot <- y[1]*(rho-y[3])-y[2]
	     zdot <- y[1]*y[2]-beta*y[3]
	     return(list(c(xdot, ydot, zdot)))
	})
}

parms <- c(sigma= 10, beta= 8/3, rho= 28)
# yini <- runif(3)
times <- seq(from=0, by= 1e-3, length.out= 1e4)
cat("running simulation\n")
system.time(
	path <- ode(func= lorenz, times= times, y=yini, parms=parms, method="rk4")
)
path <- data.frame(path)
scatterplot3d(path$X1, path$X2, path$X3, type="p", color=rainbow(n=nrow(path)), pch=".", lwd=3)

if (GIF)
{
	png("/tmp/lorenz%03d.png")
	cat("plotting to PNGs\n")
	system.time(
	for (i in round(seq(from=1, to=length(path$X1), length.out=400)))
	{
		s3d <- scatterplot3d(
				     path$X1[1:i],
				     path$X2[1:i], 
				     path$X3[1:i],
				     type="p", color=rainbow(n=nrow(path))[1:i],
				     pch= ".",
				     xlim=range(path$X1), 
				     ylim=range(path$X2), 
				     zlim=range(path$X3), 
				     )
		s3d$points3d(path$X1[i], path$X2[i], path$X3[i])
	}
	)
	dev.off()
	cat("converting PNG into GIF\n")
	system.time(system(" convert -delay 10 /tmp/lorenz*.png ./lorenz.gif "))
}
