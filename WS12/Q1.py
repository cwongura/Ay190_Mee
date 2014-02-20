# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 09:39:44 2014

For Ay190 - Computational Astrophysics - WS12 Q1

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *

data = np.loadtxt('presupernova.dat')
N = np.shape(data)[1]

for i in range(N):
    print "Column %d" % (i)
    print data[:,i]

col1 = data[:,0] # Just index
mass = data[:,1] # It is the enclosed mass.
rad = data[:,2] # Radius
temp = data[:,3] # Temperature
rho = data[:,4] # Density
rad_v = data[:,5] # Radial velocity
Ye = data[:,6] # fraction of electron per baryon

# Determine the index when the radius is going over 10^9
for i in range(len(rad)):
    n = i
    if rad[i] > 10**9:
        break

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

ax = pl.gca()
p1, = pl.plot(rad[0:n], temp[0:n], 'b', linewidth=2)
p2, = pl.plot(rad[0:n], rho[0:n], 'r', linewidth=2)

ax.set_xlabel('Radius')

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.legend((p1,p2),("Col3","Col4"),loc=(0.1,0.1),frameon=False)

pl.show()