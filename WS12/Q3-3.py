# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 11:46:33 2014

For Ay190 - Computational Astrophysics - WS12 Q2

@author: MeenoiKung
"""

import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_12 as Lib12

ggrav = 6.67e-8

data = np.loadtxt('presupernova.dat')

# We know from the previous question which column represents which quality
rad_old = data[:,2]
rho_old = data[:,4]
mass_old = data[:,1]

# Determine the index when the radius is going over 10^9
for i in range(len(rad_old)):
    n = i+1
    if rad_old[i] > 10**9:
        break
    
rad_old = rad_old[:n]
rho_old = rho_old[:n]
mass_old = mass_old[:n]

npoints = 10000
rad = np.linspace(rad_old[0],1e9, num= npoints)
dr = rad[1]-rad[0]

f = inter.interp1d(rad_old,rho_old,kind='cubic')
rho = f(rad)

# set up variables
phi = np.zeros(npoints)
z = np.zeros(npoints)
mass = np.zeros(npoints)

# set up initial values
phi[0] = 0
z[0] = 0
mass[0] = mass_old[0]

for n in range(npoints-1):
    new = Lib12.integrate_FE(rad[n],dr,phi[n],z[n],mass[n],rho[n])
    phi[n+1] = new[0]
    z[n+1] = new[1]
    mass[n+1] = new[2]

# phi at r_surface should be equal to
final = -ggrav*mass[-1]/rad[-1]

phi += final-phi[-1]

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(rad,phi,'r-',linewidth=2)

pl.xlabel('Radius')
pl.ylabel('Phi')

pl.show()