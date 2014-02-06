# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 19:57:34 2014

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

npoints = np.arange(300,4000,10)
radmax = 2.0e8
mass_con = np.zeros(len(npoints))
mass_con_RK2 = np.zeros(len(npoints))
mass_con_RK3 = np.zeros(len(npoints))
mass_con_RK4 = np.zeros(len(npoints))

for i in range(len(npoints)):
    radius = np.linspace(0, radmax, num=npoints[i])
    dr = radius[1]-radius[0]

    # set up variables
    press = np.zeros(npoints[i])
    rho   = np.zeros(npoints[i])
    mass  = np.zeros(npoints[i])
    
    press_RK2 = np.zeros(npoints[i])
    rho_RK2   = np.zeros(npoints[i])
    mass_RK2  = np.zeros(npoints[i])
    
    press_RK3 = np.zeros(npoints[i])
    rho_RK3   = np.zeros(npoints[i])
    mass_RK3  = np.zeros(npoints[i])
    
    press_RK4 = np.zeros(npoints[i])
    rho_RK4   = np.zeros(npoints[i])
    mass_RK4  = np.zeros(npoints[i])
    
    # set up central values
    rho[0]   = 1.0e10
    press[0] = polyK * rho[0]**polyG
    mass[0]  = 0.0
    
    rho_RK2[0]   = 1.0e10
    press_RK2[0] = polyK * rho[0]**polyG
    mass_RK2[0]  = 0.0
    
    rho_RK3[0]   = 1.0e10
    press_RK3[0] = polyK * rho[0]**polyG
    mass_RK3[0]  = 0.0
    
    rho_RK4[0]   = 1.0e10
    press_RK4[0] = polyK * rho[0]**polyG
    mass_RK4[0]  = 0.0
    
    # set up termination criterion
    press_min = 1.0e-10 * press[0]
    
    nsurf = 0
    for n in range(npoints[i] - 1):
    
        (press[n+1],mass[n+1]) = Lib8.tov_integrate_FE(radius[n],dr,
                                                        press[n],rho[n],mass[n])
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n

        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]

        # invert the EOS to get density
        rho[n+1] = (press[n+1]/polyK)**(1/polyG)
    
    nsurf_RK2 = 0
    for n in range(npoints[i] - 1):
    
        (press_RK2[n+1],mass_RK2[n+1]) = Lib8.tov_integrate_RK2(radius[n],dr,
                                                        press_RK2[n],rho_RK2[n],mass_RK2[n])
        # check for termination criterion
        if(press_RK2[n+1] < press_min and nsurf_RK2==0):
            nsurf_RK2 = n

        if(n+1 > nsurf_RK2 and nsurf_RK2 > 0):
            press_RK2[n+1] = press_RK2[nsurf]
            rho_RK2[n+1]   = rho_RK2[nsurf]
            mass_RK2[n+1]  = mass_RK2[nsurf]

        # invert the EOS to get density
        rho_RK2[n+1] = (press_RK2[n+1]/polyK)**(1/polyG)
        
    nsurf_RK3 = 0
    for n in range(npoints[i] - 1):
    
        (press_RK3[n+1],mass_RK3[n+1]) = Lib8.tov_integrate_RK3(radius[n],dr,
                                                        press_RK3[n],rho_RK3[n],mass_RK3[n])
        # check for termination criterion
        if(press_RK3[n+1] < press_min and nsurf_RK3==0):
            nsurf_RK3 = n

        if(n+1 > nsurf_RK3 and nsurf_RK3 > 0):
            press_RK3[n+1] = press_RK3[nsurf]
            rho_RK3[n+1]   = rho_RK3[nsurf]
            mass_RK3[n+1]  = mass_RK3[nsurf]

        # invert the EOS to get density
        rho_RK3[n+1] = (press_RK3[n+1]/polyK)**(1/polyG)
        
    nsurf_RK4 = 0
    for n in range(npoints[i] - 1):
    
        (press_RK4[n+1],mass_RK4[n+1]) = Lib8.tov_integrate_RK4(radius[n],dr,
                                                        press_RK4[n],rho_RK4[n],mass_RK4[n])
        # check for termination criterion
        if(press_RK4[n+1] < press_min and nsurf_RK4==0):
            nsurf_RK4 = n

        if(n+1 > nsurf_RK4 and nsurf_RK4 > 0):
            press_RK4[n+1] = press_RK4[nsurf]
            rho_RK4[n+1]   = rho_RK4[nsurf]
            mass_RK4[n+1]  = mass_RK4[nsurf]

        # invert the EOS to get density
        rho_RK4[n+1] = (press_RK4[n+1]/polyK)**(1/polyG)
    
    mass_con[i] = mass[nsurf]/msun
    mass_con_RK2[i] = mass_RK2[nsurf]/msun
    mass_con_RK3[i] = mass_RK3[nsurf]/msun
    mass_con_RK4[i] = mass_RK4[nsurf]/msun
    
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(npoints, np.abs(mass_con - mass_con[-1]), "r", linewidth=2)
p2, = pl.plot(npoints, np.abs(mass_con_RK2 - mass_con_RK2[-1]), "b", linewidth=2)
p3, = pl.plot(npoints, np.abs(mass_con_RK3 - mass_con_RK3[-1]), "g", linewidth=2)
p4, = pl.plot(npoints, np.abs(mass_con_RK4 - mass_con_RK4[-1])*10**(-2), "m", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of Points", labelpad = 15)
pl.ylabel("Absolute Error", labelpad = -5)

# Legend
pl.legend((p1,p2,p3,p4),("FE","RK2","RK3","RK4"),loc=(0.,0.1),frameon=False)

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()