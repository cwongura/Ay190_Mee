# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 03:53:17 2014

For Ay 190 - Computational Astrophysics - WS08

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
npoints = 5000
radmax = 2.0e8 # 2000 km
radius = np.linspace(0, radmax, num=npoints)
dr = radius[1]-radius[0]

# set up variables
press = np.zeros(npoints)
rho   = np.zeros(npoints)
mass  = np.zeros(npoints)

press_RK3 = np.zeros(npoints)
rho_RK3   = np.zeros(npoints)
mass_RK3  = np.zeros(npoints)

# set up central values
rho[0]   = 1.0e10
press[0] = polyK * rho[0]**polyG
mass[0]  = 0.0

rho_RK3[0]   = 1.0e10
press_RK3[0] = polyK * rho[0]**polyG
mass_RK3[0]  = 0.0

# set up termination criterion
press_min = 1.0e-10 * press[0]

nsurf = 0
for n in range(npoints-1):
    
    (press[n+1],mass[n+1]) = Lib8.tov_integrate_RK2(radius[n],
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

###### Uncomment the next section to calculate using RK3
#nsurf_RK3 = 0
#for n in range(npoints-1):
    
#    (press_RK3[n+1],mass_RK3[n+1]) = Lib8.tov_integrate_RK3(radius[n],
#                                              dr,
#                                              press_RK3[n],
#                                              rho_RK3[n],mass_RK3[n])
    # check for termination criterion
#    if(press_RK3[n+1] < press_min and nsurf_RK3==0):
#        nsurf_RK3 = n

#    if(n+1 > nsurf_RK3 and nsurf_RK3 > 0):
#        press_RK3[n+1] = press_RK3[nsurf]
#        rho_RK3[n+1]   = rho_RK3[nsurf]
#        mass_RK3[n+1]  = mass_RK3[nsurf]

    # invert the EOS to get density
#    rho_RK3[n+1] = (press_RK3[n+1]/polyK)**(1/polyG)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

ax1 = pl.gca()
p1, = ax1.plot(radius, mass/mass[-1], 'b', linewidth=2)
p2, = ax1.plot(radius, press/press[0], 'r', linewidth=2)
ax1.set_xlabel('Radius')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Mass or Pressure')

ax2 = ax1.twinx()
p3, = ax2.plot(radius, rho, 'g', linewidth=2)
for t1 in ax2.get_yticklabels():
    t1.set_color('g')
ax2.set_ylabel('Density', color='g', labelpad = -15)
ax2.set_yscale('Log')

pl.legend((p1,p2,p3), ("Mass/FinalMass", "Press/CentPress", "Density"),
          loc=(0.45,0.7),frameon=False)

pl.show()