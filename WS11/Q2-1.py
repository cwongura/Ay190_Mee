# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 13:12:12 2014

For Ay190 - Computational Astrophysics - WS11 Q2-1

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as mpl

def apply_bcs(x,y):
    # apply boundary conditions
    # to make sure that the derivatives are zero at both end
    y[1] = y[0]
    y[-1] = y[-2]
    return y


def analytic(x,v,t,x0,sigma):
    return np.exp(-(x-x0-v*t)**2/(2*sigma**2))

def upwind(yold,v,dt,dx):
    du = np.append(0,yold[1:]-yold[:-1])
    ynew = yold - (v*dt/dx)*du
    return ynew

# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = np.arange(0,100,0.1)
# parameters
dx = x[1]-x[0]
v = 0.1

n = len(x)
y = np.zeros(n)

# for initial data
sigma = np.sqrt(15.0)
x0 = 30.0

ntmax = 700
cfl = np.array([0.,0.5,1.,1.5])
err = np.zeros((len(cfl),ntmax))

for m in range(len(cfl)):
    t = 0.
    dt = cfl[m]*dx/np.abs(v)
    
    #set up initial conditions
    y = analytic(x,v,t,x0,sigma)
    
    yold = y
    for it in range(ntmax):
        t = (it+1)*dt
        
        # save previous and previous previous data
        yold = y
        # get new data; ideally just call a function
        y = upwind(yold, v, dt, dx)
        
        # after update, apply boundary conditions
        # apply_bcs(x,y) 
        y = apply_bcs(x,y)

        # get analytic result for time t
        yana = analytic(x,v,t,x0,sigma)
        
        err[m,it] = np.max(np.abs(yana-y))


p1, = mpl.plot(range(ntmax),err[0,:],'r-')
p2, = mpl.plot(range(ntmax),err[1,:],'b-')
p3, = mpl.plot(range(ntmax),err[2,:],'g-')
p4, = mpl.plot(range(ntmax),err[3,:],'m-')


mpl.legend((p1,p2,p3,p4),('alpha=0.','alpha=0.5','alpha=1.',
           'alpha=1.5'),loc=(0.,0.75),frameon=False)
           
mpl.xlabel('Step')
mpl.ylabel('Error')

mpl.xscale('log')
mpl.yscale('log')

mpl.show()