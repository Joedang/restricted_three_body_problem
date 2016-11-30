saveVideo({
	ani.options(interval = 0.05, nmax = 300)
	i=0
	for (i in 1:300)
		plot(i)
}, video.name = "BM.mp4", other.opts = "-pix_fmt yuv420p -b 300k")

