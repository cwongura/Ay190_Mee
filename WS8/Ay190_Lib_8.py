#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl
from plot_defaults import *


# global constants
ggrav = 6.67e-8
msun  = 1.99e33

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG


#######################################
# function definitions
def tov_RHS(rad,p,rho,m):
    
    # RHS function
    
    rhs = np.zeros(2)
    if(rad > 1.0e-10):
        rhs[0] = -ggrav*m*rho/rad**2
        rhs[1] = 4*np.pi*rho*rad**2
    else:
        rhs[0] = 0.0
        rhs[1] = 0.0

    return rhs

def tov_integrate_FE(rad,dr,p,rho,m):

    # Forward-Euler Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m

    # forward Euler integrator
    new = old + dr*tov_RHS(rad, p, rho, m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew)

def tov_integrate_RK2(rad,dr,p,rho,m):

    # RK2 Integrator

    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    
    # Calculate constants     
    
    k1 = dr*tov_RHS(rad, p, rho, m)
    rho1 = ((p+0.5*k1[0])/polyK)**(1/polyG)
    k2 = dr*tov_RHS(rad+0.5*dr, p+0.5*k1[0], rho1, m+0.5*k1[1])
    
    # Perform the RK2 method
    new = old + k2
    
    return (new[0], new[1])

def tov_integrate_RK3(rad, dr, p, rho, m):
    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    
    # Calculate constants
    k1 = dr*tov_RHS(rad, p, rho, m)
    rho1 = ((p+0.5*k1[0])/polyK)**(1/polyG)
    k2 = dr*tov_RHS(rad+0.5*dr, p+0.5*k1[0], rho1, m+0.5*k1[1])
    rho2 = ((p-k1[0]+2*k2[0])/polyK)**(1/polyG)
    k3 = dr*tov_RHS(rad+dr,p-k1[0]+2*k2[0],rho2,m-k1[1]+2*k2[1])
    
    # Perform the RK3 method
    new = old + (1./6)*(k1+4*k2+k3)
    
    return (new[0], new[1])

def tov_integrate_RK4(rad, dr, p, rho, m):
    new = np.zeros(2)
    old = np.zeros(2)
    old[0] = p
    old[1] = m
    
    # Calculate constants
    k1 = dr*tov_RHS(rad, p, rho, m)
    rho1 = ((p+0.5*k1[0])/polyK)**(1/polyG)
    k2 = dr*tov_RHS(rad+0.5*dr, p+0.5*k1[0], rho1, m+0.5*k1[1])
    rho2 = ((p+0.5*k2[0])/polyK)**(1/polyG)
    k3 = dr*tov_RHS(rad+0.5*dr, p+0.5*k2[0], rho2, m+0.5*k2[1])
    rho3 = ((p+k3[0])/polyK)**(1/polyG)
    k4 = dr*tov_RHS(rad+dr, p+k3[0], rho3, m+k3[1])
    
    # Perform the RK4 method
    new = old + (1./6)*(k1+2*k2+2*k3+k4)
    
    return (new[0], new[1])
    