# -*- coding: utf-8 -*-
"""
For Ay 190 - Computational Astrophysics - Final Project

This python code explores the three-body system: Sun-Earth-Moon. It produces
two figures.
The first figure shows the distance between Earth and Moon over time. We can
see how the perigee and apogee varies over time.
The second figure measures eccentricity of the moon for each orbit.

@author: MeenoiKung
"""

import numpy as np
import Integrator as Int
import Project_Lib as Lib
import matplotlib.pyplot as pl
from plot_defaults import *

# global constants
ggrav = 6.67e-8
msun = 1.99e33
mearth = 5.97e27
mmoon = 7.342e25
au = 1.5e13
dmoon = 3.844e10
eearth = 0.0167
emoon = 0.0549
seconds_per_year = 24.*3600*365

# system parameters
Nsteps = 30000
t0 = 0
t1 = 2 * seconds_per_year
dt = (t1-t0)/Nsteps
t = np.linspace(t0, t1, num = Nsteps)

# setup initial parameters for the system
mass = np.array([msun, mearth, mmoon])
x = np.zeros(3)
y = np.zeros(3)
z = np.zeros(3)
vx = np.zeros(3)
vy = np.zeros(3)
vz = np.zeros(3)

# Put both Earth and moon at perigee
x[0] = 0.*au # Sun
x[1] = 1.*au*(1-eearth) # Earth
x[2] = x[1] + dmoon*(1-emoon)

# Select initial velocity that makes Earth orbits the sun counterclockwise
# and the moon orbits the earth counterclockwise
# The magnitude of velocity would put both Earth and moon in eccentric orbit
# around the sun and earth, respectively, in the two-body system.
vy[0] = 0.
vy[1] = np.sqrt((ggrav*msun/au)*(1+eearth)/(1-eearth))
vy[2] = np.sqrt((ggrav*mearth/dmoon)*(1+emoon)/(1-emoon)) + vy[1]

u = np.array((x,y,z,vx,vy,vz)).transpose()

# An array collecting distance between Earth-Moon.
dist = np.zeros(Nsteps)
dist[0] = Lib.CalDistEM(u)

for it in range(1, Nsteps):
    print it
    time = t0 + it*dt
    u = Int.NbodyRK4(u, mass, time, dt)
    
    dist[it] = Lib.CalDistEM(u)

p1, = pl.plot(t/seconds_per_year, dist, 'b-')
pl.xlabel("Time (Years)")
pl.ylabel("Distance Earth-Moon")

# Use helper function to calculate eccentricity
ecc = Lib.CalEccMoon(dist, t)
pl.figure()
p2, = pl.plot(np.arange(len(ecc)), ecc, 'r-')
pl.xlabel("Cycle")
pl.ylabel("Eccentricity")

pl.show()