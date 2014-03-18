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
    
def CalDistEM(u):
    '''
    Calculate distance between Earth and Moon
    '''
    xearth = u[1,0]
    xmoon = u[2,0]
    yearth = u[1,1]
    ymoon = u[2,1]
    zearth = u[1,2]
    zmoon = u[2,2]
    
    return np.sqrt((xearth-xmoon)**2+(yearth-ymoon)**2+(zearth-zmoon)**2)

def CalEccMoon(dist, time):
    '''
    Calculate eccentricity based on apogee and perigee of each cycle.
    Assuming the orbital period of the moon is 27 days.
    '''
    period = 27*24*3600 # 27 days in second
    t0 = time[0]
    t1 = t0 + period
    ini = 0
    
    ecc = np.array([])
    for i in range(len(time)):
        if time[i] > t1:
            r_a = np.max(dist[ini:i])
            r_p = np.min(dist[ini:i])
            ecc = np.append(ecc, (r_a-r_p)/(r_a+r_p))
            
            ini = i
            t1 += period
        if t1 > time[-1]:
            break
    return ecc