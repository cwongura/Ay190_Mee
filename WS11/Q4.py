# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:05:04 2014

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
    
def LaxFried(yold, v, dt, dx):
    plus = yold[2:]+yold[:-2]
    minus = yold[2:]-yold[:-2]
    ynew = 0.5*plus - (0.5*v*dt/dx)*minus
    ynew = np.append(yold[0],ynew)
    ynew = np.append(ynew,yold[-1])
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
sigma = np.sqrt(15.0)
x0 = 30.0

#set up initial conditions
y_uw = analytic(x,v,t,x0,sigma)
y_lf = analytic(x,v,t,x0,sigma)

# evolve (and show evolution)
mpl.ion()
mpl.figure()
mpl.plot(x,y,'b-') # numerical data
mpl.plot(x,analytic(x,v,t,x0,sigma),'r-') # analytic data
mpl.show()


yold_uw = y_uw
yold_lf = y_lf
ntmax = 1300
for it in range(ntmax):
    t = (it+1)*dt
    # save previous and previous previous data
    yold_uw = y_uw
    yold_lf = y_lf

    # get new data; ideally just call a function
    y_uw = upwind(yold_uw, v, dt, dx)
    y_lf = LaxFried(yold_lf, v, dt, dx)

    # after update, apply boundary conditions
    # apply_bcs(x,y) 
    y_uw = apply_bcs(x,y_uw)
    y_lf = apply_bcs(x,y_lf)

    # get analytic result for time t
    yana = analytic(x,v,t,x0,sigma)
    # compute error estimage
    err = 0

    err = np.max(np.abs(y_uw - y_lf))
    print "it = ",it
    print "err =", err
    mpl.clf()

    # plot results
    p1, = mpl.plot(x,yana,'r-')
    p2, = mpl.plot(x,y_uw, 'b-')
    p3, = mpl.plot(x,y_lf, 'g-')
    mpl.text(0.5,0.97,'t=%f' % t)
    mpl.text(0.5,0.9,err)
    mpl.legend((p1,p2,p3),('Analytical','Upwind','Lax-Friedrich'),
               loc=(0.,0.),frameon=False)
    mpl.draw()

mpl.show()