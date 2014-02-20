# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 06:05:25 2014

For Ay 190 - Computational Astrophysics - WS 12 - Q3-1

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_12 as Lib12
from plot_defaults import *


ggrav = 6.67e-8

rmin = 0.
rmax = 1e9
npoints = 1000

# set up grid
radius = np.linspace(rmin, rmax, num = npoints)
dr = radius[1]-radius[0]

# set up variables
phi = np.zeros(npoints)
z = np.zeros(npoints)
mass = np.zeros(npoints)

# set up initial values
phi[0] = 0
z[0] = 0
mass[0] = 0

# set the density to be constant 1e10
rho = np.zeros(npoints) + 1e10

for n in range(npoints-1):
    new = Lib12.integrate_FE(radius[n],dr,phi[n],z[n],mass[n],rho[n])
    phi[n+1] = new[0]
    z[n+1] = new[1]
    mass[n+1] = new[2]

# phi at r_surface should be equal to
final = -ggrav*mass[-1]/radius[-1]

phi += final-phi[-1]

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

analytic = (2./3)*np.pi*ggrav*rho[0]*(radius**2-3*rmax**2)

p1, = pl.plot(radius, phi, 'r-', linewidth=2)
p2, = pl.plot(radius, analytic,'b-',linewidth=2)

# label the axes
pl.xlabel('Radius', labelpad = 15)
pl.ylabel('phi', labelpad=-5)

pl.legend((p1,p2),('Numerical','Analytic'),loc=(0.,0.8),frameon=False)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p3, = pl.plot(radius, np.abs(phi-analytic),'g-',linewidth=2)

pl.xscale('log')
pl.yscale('log')

pl.show()