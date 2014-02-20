# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 06:48:16 2014

For Ay 190 - Computational Astrophysics - WS 12 - Q3-2

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_12 as Lib12
from plot_defaults import *

ggrav = 6.67e-8

rmin = 0.
rmax = 1e9
npoints = np.array([100,1000,10000,100000])
err = np.zeros(len(npoints))

for i in range(len(npoints)):
    # set up grid
    radius = np.linspace(rmin, rmax, num = npoints[i])
    dr = radius[1]-radius[0]
    
    # set up variables
    phi = np.zeros(npoints[i])
    z = np.zeros(npoints[i])
    mass = np.zeros(npoints[i])
    
    # set up initial values
    phi[0] = 0
    z[0] = 0
    mass[0] = 0
    
    # set the density to be constant 1e10
    rho = np.zeros(npoints[i]) + 1e10
    
    for n in range(npoints[i]-1):
        new = Lib12.integrate_FE(radius[n],dr,phi[n],z[n],mass[n],rho[n])
        phi[n+1] = new[0]
        z[n+1] = new[1]
        mass[n+1] = new[2]
    
    # phi at r_surface should be equal to
    final = -ggrav*mass[-1]/radius[-1]

    phi += final-phi[-1]
    analytic = (2./3)*np.pi*ggrav*rho[0]*(radius**2-3*rmax**2)
    
    err[i] = np.abs(phi[-1]-analytic[-1])

p1, = pl.plot(npoints,err,'r-',linewidth=2)

pl.xlabel('Numer of points')
pl.ylabel('Absolute error')

pl.xscale('log')
pl.yscale('log')

pl.show()