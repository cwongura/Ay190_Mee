# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 19:21:35 2014

For Ay 190 - Computation Astrophysics - WS08

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_8 as Lib8

# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

# set up grid
npoints = 1000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0, radmax, num=npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = Lib8.tov_integrate_FE(radius[n],
                                              dr,
                                              press[n],
                                              rho[n],mass[n])
    # check for termination criterion
    if(press[n+1] < press_min and nsurf==0):
        nsurf = n

    if(n+1 > nsurf and nsurf > 0):
        press[n+1] = press[nsurf]
        rho[n+1]   = rho[nsurf]
        mass[n+1]  = mass[nsurf]

    # invert the EOS to get density
    rho[n+1] = (press[n+1]/polyK)**(1/polyG)


print radius[nsurf]/1.0e5
print mass[nsurf]/msun