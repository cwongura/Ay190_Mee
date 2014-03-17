# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 04:17:24 2014

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
au = 1.5e13
seconds_per_year = 24.*3600*365

# system parameters
Nsteps = 3000
t0 = 0
t1 = 80 * seconds_per_year
dt = (t1-t0)/Nsteps
t = np.linspace(t0, t1, num = Nsteps)

# setup initial parameters for the system
mass = np.array([msun, mearth])
x = np.zeros(2)
y = np.zeros(2)
z = np.zeros(2)
vx = np.zeros(2)
vy = np.zeros(2)
vz = np.zeros(2)

x[0] = -(mearth/(msun+mearth))*au
x[1] = (msun/(msun+mearth))*au

# Set up velocity to ensure that it's a circular orbit
vy[0] = -np.sqrt(ggrav*mearth*np.abs(x[0])/au**2)
vy[1] = np.sqrt(ggrav*msun*x[1]/au**2)

uRK4 = np.array((x,y,z,vx,vy,vz)).transpose()
uLF = np.array((x,y,z,vx,vy,vz)).transpose()

energy_RK4 = np.zeros(Nsteps)
energy_LF = np.zeros(Nsteps)

energy_RK4[0] = Lib.CalTotalE(uRK4, mass)
energy_LF[0] = Lib.CalTotalE(uLF, mass)

for it in range(1, Nsteps):
    print it
    time = t0 + it*dt
    uRK4 = Int.NbodyRK4(uRK4, mass, time, dt)
    uLF = Int.NbodyLeapFrog(uLF, mass, dt)

    energy_RK4[it] = Lib.CalTotalE(uRK4, mass)
    energy_LF[it] = Lib.CalTotalE(uLF, mass)

energy_RK4 = np.abs((energy_RK4-energy_RK4[0])/energy_RK4[0])
energy_LF = np.abs((energy_LF-energy_LF[0])/energy_LF[0])

p1, = pl.plot(t/seconds_per_year, energy_RK4, 'r-')
p2, = pl.plot(t/seconds_per_year, energy_LF, 'b-')
pl.xlabel("Time(year)")
pl.ylabel("dE/E")
pl.yscale('log')
pl.legend((p1,p2),("RK4","LF"),loc=(0.,0.),frameon=False)
pl.show()