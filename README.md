# pipe-conductance
Computes conductance of pipes of arbitrary but uniform cross section in the molecular flow regime. This regime is also called the rarified gas regime. It means that collisions between particles are more rare than collisions between particles and the pipe. More technically, the mean free path of a particle in the gas should be longer than the diameter or other relevant length scale of the pipe. 

In practice, this regime is relevant in vacuum chambers at room temperature when the pressure is below about a militorr, a million times less pressure than atmospheric pressure, or a tenth of a Pascal. The size of the chamber and the mass of the gas flowing in it will effect this.

The code is implemented in MATLAB, and thus can only be run with MATLAB and a MATLAB license. It could be ported to another language rather easily, with the exception of the handling of arbitrary pipe cross sections specified as a black and white image array. This would require an image processing package for extracting the boundary from the image, but may be worth the effort at some point.

