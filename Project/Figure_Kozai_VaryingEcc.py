# -*- coding: utf-8 -*-
"""
For Ay 190 - Computational Astrophysics - Final Project

Setting the moon with initial inclination of 45 degree to the Earth-Sun
orbital plane. The code produces figure that observe moon eccentricity.

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
seconds_per_year = 24.*3600*365

theta = 45 * (np.pi/180.)

# system parameters
Nsteps = 300000
t0 = 0
t1 = 200 * seconds_per_year
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
x[1] = 1.*au # Earth
x[2] = x[1] + dmoon*np.cos(theta) # Put moon at theta inclination
z[2] = dmoon*np.sin(theta)

vy[0] = 0.
vy[1] = np.sqrt(ggrav*msun/au)
vy[2] = np.sqrt(ggrav*mearth/dmoon) + vy[1]

u = np.array((x,y,z,vx,vy,vz)).transpose()

dist = np.zeros(Nsteps)
dist[0] = Lib.CalDistEM(u)

for it in range(1, Nsteps):
    print it
    time = t0 + it*dt
    u = Int.NbodyRK4(u, mass, time, dt)

    dist[it] = Lib.CalDistEM(u)

# Use helper function to calculate eccentricity
ecc = Lib.CalEccMoon(dist, t)

p2, = pl.plot(np.arange(len(ecc)), ecc, 'r-')
pl.xlabel("Cycle")
pl.ylabel("Eccentricity")

pl.show()