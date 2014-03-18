# -*- coding: utf-8 -*-
"""
For Ay 190 - Computational Astrophysics - Final Project

This code creates two figures
The first figure shows error evolution of RK4 with different time step.
The second figure shows convergence of error.

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
Nsteps = np.array([1000,2000,3000,4000,5000,6000,50000])
t0 = 0
t1 = 80 * seconds_per_year
for m in range(len(Nsteps)):
    dt = (t1-t0)/Nsteps[m]
    t = np.linspace(t0, t1, num = Nsteps[m])

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

    energy_RK4 = np.zeros(Nsteps[m])

    energy_RK4[0] = Lib.CalTotalE(uRK4, mass)
    
    for it in range(1, Nsteps[m]):
        print it
        time = t0 + it*dt
        uRK4 = Int.NbodyRK4(uRK4, mass, time, dt)

        energy_RK4[it] = Lib.CalTotalE(uRK4, mass)

    energy_RK4 = np.abs((energy_RK4-energy_RK4[0])/energy_RK4[0])
    
    if m == 0:
        err_1000 = energy_RK4
        dt_1000 = dt
        t_1000 = t
    if m == 1:
        err_2000 = energy_RK4
        dt_2000 = dt
        t_2000 = t
    if m == 2:
        err_3000 = energy_RK4
        dt_3000 = dt
        t_3000 = t
    if m == 3:
        err_4000 = energy_RK4
        dt_4000 = dt
        t_4000 = t
    if m == 4:
        err_5000 = energy_RK4
        dt_5000 = dt
        t_5000 = t
    if m == 5:
        err_6000 = energy_RK4
        dt_6000 = dt
        t_6000 = t
    if m == 6:
        err_small = energy_RK4
        dt_small = dt
        t_small = t

p1, = pl.plot(t_1000/seconds_per_year, err_1000, 'r-')
p2, = pl.plot(t_2000/seconds_per_year, err_2000, 'b-')
p3, = pl.plot(t_3000/seconds_per_year, err_3000, 'g-')
p4, = pl.plot(t_4000/seconds_per_year, err_4000, 'm-')
p5, = pl.plot(t_5000/seconds_per_year, err_5000, 'k-')
p6, = pl.plot(t_6000/seconds_per_year, err_6000, 'c-')
p7, = pl.plot(t_small/seconds_per_year, err_small, 'r--')

pl.xlabel("Time(year)")
pl.ylabel("dE/E")
pl.yscale('log')
pl.legend((p1,p2,p3,p4,p5,p6,p7),("dt="+str(dt_1000),"dt="+str(dt_2000),
"dt="+str(dt_3000),"dt="+str(dt_4000),"dt="+str(dt_5000),"dt="+str(dt_6000),
"dt="+str(dt_small)),
loc=(0.77,0.),frameon=True)

pl.figure()
p8, = pl.plot(np.array([dt_1000,dt_2000,dt_3000,dt_4000,dt_5000,dt_6000,
                        dt_small]),np.array([np.max(err_1000),np.max(err_2000),
np.max(err_3000),np.max(err_4000),np.max(err_5000),np.max(err_6000),
np.max(err_small)]),'b-')
pl.xlabel("Time step")
pl.ylabel("Maximum absolute error")
pl.xscale('log')
pl.yscale('log')

pl.show()