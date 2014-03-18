# -*- coding: utf-8 -*-
"""
For Ay 190 - Computational Astrophysics - Final Project

Creates figure that shows error in energy of circular and eccentric orbit, with
the same size of time step.

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
e = 0.5

# system parameters
Nsteps = 2500
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

# Set up for circular orbit
x[1] = au
vy[1] = np.sqrt(ggrav*msun/au)

u_cir = np.array((x,y,z,vx,vy,vz)).transpose()

# Set up for eccentric orbit
x[1] = au*(1-e)
vy[1] = np.sqrt((ggrav*msun/au)*(1+e)/(1-e))

u_ecc = np.array((x,y,z,vx,vy,vz)).transpose()

energy_cir = np.zeros(Nsteps)
energy_ecc = np.zeros(Nsteps)

energy_cir[0] = Lib.CalTotalE(u_cir, mass)
energy_ecc[0] = Lib.CalTotalE(u_ecc, mass)

for it in range(1, Nsteps):
    print it
    time = t0 + it*dt
    u_cir = Int.NbodyRK4(u_cir, mass, time, dt)
    u_ecc = Int.NbodyRK4(u_ecc, mass, time, dt)

    energy_cir[it] = Lib.CalTotalE(u_cir, mass)
    energy_ecc[it] = Lib.CalTotalE(u_ecc, mass)

energy_cir = np.abs((energy_cir-energy_cir[0])/energy_cir[0])
energy_ecc = np.abs((energy_ecc-energy_ecc[0])/energy_ecc[0])

p1, = pl.plot(t/seconds_per_year, energy_cir, 'r-')
p2, = pl.plot(t/seconds_per_year, energy_ecc, 'b-')
pl.xlabel("Time(year)")
pl.ylabel("dE/E")
pl.yscale('log')
pl.legend((p1,p2),("Circular","Eccentric"),loc=(0.,0.),frameon=False)
pl.show()