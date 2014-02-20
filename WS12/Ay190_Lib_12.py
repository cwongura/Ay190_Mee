# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 06:13:58 2014

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl

ggrav = 6.67e-8

def RHS(rad, phi, z, mass, rho):
    if rad==0.:
        rad = 1e-16
    
    dphi = z
    dz = 4*np.pi*ggrav*rho-2*z/rad
    dm = 4*np.pi*(rad**2)*rho
    
    return np.array([dphi,dz,dm])

def integrate_FE(rad,dr,phi,z,mass,rho):
    old = np.zeros(3)
    new = np.zeros(3)
    
    old[0] = phi
    old[1] = z
    old[2] = mass
    
    new = old + dr*RHS(rad, phi, z, mass, rho)
    
    return new