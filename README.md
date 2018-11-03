# restricted_three_body_problem
Some tools for exploring the restricted 3 body problem with R.

[If you're looking for pretty pictures, click here!](interestingOrbits/confinedTransfers.mp4)  
_Shown: a satellite that's bound by the energy-like potential in the synodic reference frame._

I originally did this as part of a report on the 3 body problem for a PhD-level course in classical mechanics. In truth, it was largely unnecessary for the repot and just an excuse to write/play with a simulation.
The bulk of the work was on figuring out how to set up certain scenarios and automating a whole bunch of visualizations, so I could tell what kind of orbit I was actually looking at.

I focused heavily on the synodic coordinate system (rotational reference frame locked to the two larger masses), so some of the plots labeled "kinetic energy" are actually one term of the energy-like conserved quantity in this coordinate system. 

This is also the first time I've experimented with having a simulation that automatically journals its settings and results, which is incredibly useful for when I want to go back to a scenario that I was experimenting with but didn't save.

### How to Use
Tweak the input parameters in orbit.R to look at different scenarios. 
You can load the environments from `interestingOrbits/` to view some interesting scenarios. 
Any time `orbit.R` is run, the environment is automatically saved to `interestingOrbits/archive/`, which is useful if you forget how you set up a particular scenario. 
Tweak the function calls in orbitPlot.R to change the output plots. 
Call `GIF <- T` before running orbitPlot.R to export a GIF of the orbital motion. 
This feature actually creates PDFs for each frame and then converts them to JPG and builds those into a GIF using imagemagick (via orbitGIF.sh).
I did this because I did all the development on my Chromebook, on which I couldn't get R to work with `png()` or `jpg()`, which are required to make GIFs in the usual libraries.
If you have an installation of R where those functions work, you can just replace my GIF-making method with like one function call from one of the existing packages that do that.

### depends
* R, obviously
  * `deSolve` library
* If you're using my method for creating GIFs:
  * imagemagick
  * bash
  * Unix-like system with `/tmp/`

