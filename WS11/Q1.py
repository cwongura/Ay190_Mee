# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:04:51 2014

For Ay 190 - Computational Astrophysics - WS11 Q1

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as mpl

def analytic(x,v,t,x0,sigma):
    return np.exp(-(x-x0-v*t)**2/(2*sigma**2))
    
# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.arange(0,100,0.1)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)

cfl = 1.0
dt = cfl*dx/np.abs(v)
t = 0.0

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.text(0.5,0.97,'t=%f s' % t)
mpl.plot(x,analytic(x,v,t,x0,sigma),'r-') # analytic data
mpl.show()

ntmax = 700
for it in range(ntmax):
    t = (it+1)*dt
    print "it = ",it
    mpl.clf()
    mpl.plot(x,analytic(x,v,t,x0,sigma),'r-')
    mpl.text(0.5,.97,'t=%f' % t)
    mpl.pause(0.0001)
    
mpl.show()
    