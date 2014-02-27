# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 00:21:25 2014

For Ay190 - Computational Astrophysics - WS 14

@author: MeenoiKung
"""

import sys,math
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp

def apply_bcs(x,y):
    # apply boundary conditions
    # to make sure that the derivatives are zero at both end
    y[1] = y[0]
    y[-1] = y[-2]
    return y


# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
L = 100

x = np.arange(0,L,0.1)
# parameters
dx = x[1]-x[0]

n = len(x)
y = np.zeros(n)

dt = 0.5
t = 0.0

#set up initial conditions
y = (1./8)*np.sin(2.*np.pi*x/L)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'b-') # numerical data
mpl.show()

yold = y
ntmax = 400
for it in range(ntmax):
    t = (it+1)*dt
    # save previous and previous previous data
    yold = y
    ynew = np.zeros(n)
    
    for i in range(1,n-1):        
        # get new data; ideally just call a function
        du = np.append(0,yold[1:]-yold[:-1])
        if yold[i] > 0:
            ynew[i] = yold[i] - (yold[i]*dt/dx)*(yold[i]-yold[i-1])
        else:
            ynew[i] = yold[i] - (yold[i]*dt/dx)*(yold[i+1]-yold[i])
    
    y = ynew
    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y = apply_bcs(x,y)

    print "it = ",it
    print "time =", t
    mpl.clf()
    # plot numerical result
    mpl.plot(x,y,'b-')
    mpl.draw()

mpl.show()