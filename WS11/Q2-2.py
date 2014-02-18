# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 13:41:25 2014

Ay190 - Computational Astrophysics - WS11 - Q2-2

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
cfl = 0.5
dt = cfl*dx/np.abs(v)
t = 0.0

# for initial data
sigma = np.sqrt(15.0)/5.
x0 = 30.0

#set up initial conditions
y = analytic(x,v,t,x0,sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'b-') # numerical data
mpl.plot(x,analytic(x,v,t,x0,sigma),'r-') # analytic data
mpl.show()

yold2 = y
yold = y
ntmax = 700
for it in range(ntmax):
    t = (it+1)*dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    y = upwind(yold, v, dt, dx)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y = apply_bcs(x,y)

    # get analytic result for time t
    yana = analytic(x,v,t,x0,sigma)
    # compute error estimage
    err = 0

    err = np.max(np.abs(yana-y))
    print "it = ",it
    print "err =", err
    mpl.clf()
    # plot numerical result
    mpl.plot(x,y,'b-')
    # plot analytic results
    mpl.plot(x,yana,'r-')
    mpl.text(0.5,0.97,'t=%f' % t)
    mpl.text(0.5,0.9,err)
    mpl.draw()


mpl.show()