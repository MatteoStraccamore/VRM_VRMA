# Velocity Ring Model (VRM)
This code implements a VRM for analyzing velocity fields of external galaxies. The code can be used to fit a variable number of arcs (N_a = 1, 2, 4, ...) to a given galaxy and obtain two-dimensional maps of tangential and radial velocities (v_t(R, theta) and v_r(R, theta)).
## Installation
This code requires the following Python libraries: numpy, matplotlib, and scipy. These can be installed via pip:
```
pip install numpy
pip install matplotlib
pip install scipy
```
## Usage
To use the VRM, instantiate the VRM class with the following arguments:
1. Vreal: a numpy array containing the galaxy's velocity field data, with columns in the order xp, yp, and Vlos (line-of-sight velocity).
2. Iinput: the inclination angle of the galaxy in degrees.
3. phi0input: the position angle of the galaxy in degrees.
4. number_rings: the number of rings to divide the galaxy into.
5. number_arch: the number of arcs to fit to the galaxy.
6. plot: a boolean flag indicating whether to plot the generated velocity maps (default is False).
7. save: a boolean flag indicating whether to save the generated velocity maps (default is False).

Once instantiated, the VRM object can be used to generate velocity maps by calling the _matrix() method, which returns a numpy array containing the tangential and radial velocity values for each ring and arc.
The _calculation_VLOS(predictions, file) method can be used to generate a plot of the real and generated velocity maps, as well as to save the generated velocity map data to a file.
## Example
```
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
from vrm import VRM

# Load velocity field data
data = np.loadtxt('galaxy_data.txt')
vrm = VRM(data, Iinput=45, phi0input=30, number_rings=4, number_arch=1, plot=True, save=True)

# Calculate tangential and radial velocities
predictions = vrm._matrix()

# Generate velocity map plot and save data
vrm._calculation_VLOS(predictions, file='galaxy')
```
