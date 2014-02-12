#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# set up grid
xmin = 0.0
xmax = 1.0
npoints = np.array([1e1,1e2,1e3,1e4,1e5])
err_FE = np.zeros(5)
error_RK2 = np.zeros(5)

# boundary values
A = 0. # inner boundary
B = 0.1 # outer boundary

def calc_rhs(u,xx):
    # rhs routine
    # rhs[0] is rhs for y
    # rhs[1] is rhs for u
    rhs = np.zeros(2)
    rhs[0] = u
    rhs[1] = 12.*xx - 4.

    return rhs

def integrate_FE(z,x,n):
    # forward-Euler integrator
    
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    n = np.int(n)
    yy = np.zeros((n,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(n-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])

    return yy

def integrate_RK2(z,x,n):
    # RK2 integrator

    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    n = np.int(n)
    yy = np.zeros((n,2))
    
    yy[0,0] = A # boundary value A for y at x = 0
    yy[0,1] = z # guessed boundary value for y' at x = 0
    
    for i in range(n-1):
        k1 = dx*calc_rhs(yy[i,1],x[i])
        k2 = dx*calc_rhs(yy[i,1]+0.5*k1[1],x[i]+0.5*dx)
        yy[i+1,:] = yy[i,:] + k2

    return yy

for n in range(5):
    # set up grid
    x = np.linspace(xmin, xmax, num = npoints[n])
    # dx based on x[1] and x[0]
    dx = x[1] - x[0]
    
    # get initial guess for derivative
    z0 = -1100000.0
    z1 = -10000000.0
    z0_RK2 = -1100000.0
    z1_RK2 = -10000000.0
    yy0 = integrate_FE(z0,x,npoints[n])
    yy1 = integrate_FE(z1,x,npoints[n])
    yy0_RK2 = integrate_RK2(z0,x,npoints[n])
    yy1_RK2 = integrate_RK2(z1,x,npoints[n])
    phi0 = yy0[npoints[n]-1,0] - B
    phi1 = yy1[npoints[n]-1,0] - B
    phi0_RK2 = yy0_RK2[npoints[n]-1,0] - B
    phi1_RK2 = yy1_RK2[npoints[n]-1,0] - B
    dphidz = (phi1-phi0)/(z1-z0) # dphi/dz
    dphidz_RK2 = (phi1_RK2-phi0_RK2)/(z1_RK2-z0_RK2)

    i = 0
    itmax = 100
    err = 1.0e99
    err_RK2 = err
    criterion = 1.0e-12

    z0 = z1
    phi0 = phi1
    z0_RK2 = z1_RK2
    phi0_RK2 = phi1_RK2

    while (err > 1.0e-12 and i < itmax):
        z1 = z0 - phi0/dphidz # secand update
        yy = integrate_FE(z1,x,npoints[n])
        phi1 = yy[npoints[n]-1,0] - B
        dphidz = (phi1-phi0)/(z1-z0) # dphi/dz numerical
        err = np.abs(phi1) # your error measure
        z0 = z1
        phi0 = phi1
        i = i+1

        print i,z1,phi1

    i = 0
    while (err_RK2 > 1.0e-12 and i < itmax):
        z1_RK2 = z0_RK2 - phi0_RK2/dphidz_RK2 # secand update
        yy_RK2 = integrate_RK2(z1_RK2,x,npoints[n])
        phi1_RK2 = yy_RK2[npoints[n]-1,0] - B
        dphidz_RK2 = (phi1_RK2-phi0_RK2)/(z1_RK2-z0_RK2) # dphi/dz numerical
        err_RK2 = np.abs(phi1_RK2) # your error measure
        z0_RK2 = z1_RK2
        phi0_RK2 = phi1_RK2
        i = i+1

        print i,z1_RK2,phi1_RK2

    err_FE[n] = np.max(np.abs(yy[:,0]-(2.0*x**3 - 2*x**2 + 0.1*x)))
    error_RK2[n] = np.max(np.abs(yy_RK2[:,0]-(2.0*x**3 - 2*x**2 + 0.1*x)))
    
p1, = plt.plot(npoints,err_FE,"r-")
p2, = plt.plot(npoints, error_RK2, "b-")

plt.legend((p1,p2),("FE Method", "RK2 Method"),loc=(0.,0.3),frameon=False)
plt.xlabel("Number of points")
plt.ylabel("Absolute Error")
plt.xscale("Log")
plt.yscale("Log")
plt.show()