# -*- coding: utf-8 -*-
"""
For Ay 190 - Final Project - Library

@author: MeenoiKung
"""

import numpy as np

def CalUg(u, mass):
    '''
    Calculate the total gravitational potential energy of the system
    '''
    x = u[:,0]
    y = u[:,1]
    z = u[:,2]
    
    n = len(x)
    Ug = 0.
    
    for i in range(n-1):
        x_new = x[i+1:]
        y_new = y[i+1:]
        z_new = z[i+1:]
        m_new = mass[i+1:]
        
        r = np.sqrt((x[i]-x_new)**2 + (y[i]-y_new)**2 + (z[i]-z_new)**2)
        Ug -= np.sum(mass[i]*m_new/r)
    
    return Ug

def CalEk(u, mass):
    '''
    Calculate the total kinetic energy of the system
    '''
    vx = u[:,3]
    vy = u[:,4]
    vz = u[:,5]
    
    energy = 0.5*mass*(vx*vx + vy*vy + vz*vz)    
    
    return np.sum(energy)

def CalTotalE(u, mass):
    '''
    Calculate total energy of the system
    '''
    return CalUg(u, mass)+CalEk(u, mass)