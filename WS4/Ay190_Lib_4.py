# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 07:33:48 2014

Library for Ay190 Computational Astrophysics - WS 04

@author: Chatarin (Mee) Wong-u-railertkun
"""
#--------------------------------------------
import numpy as np
        
#--------------------------------------------
# Work Sheet 4
# Roots Finding
#--------------------------------------------

def Efunc(t, T, e, E):
    w = 2*np.pi/T  # angular velocity
    return E - w*t - e*np.sin(E)

def Edfunc(e, E):
    return 1 - e*np.cos(E)

def Newton(func, dfunc, x, error, t, e = 0.0167, T = 365.25635):
    n = 0
    
    if func(t, T, e, x) == 0.: # If it's already a zero, we have our answer.
        return x, func(t, T, e, x)
    while True:
        dx = func(t, T, e, x)/dfunc(e, x)
        #if dx > 0.:
        #    dx = min(dx, 10**(-2))
        #else:
        #    dx= max(dx, -10**(-2))
        x1 = x - dx        
        relerr = np.abs((func(t, T, e, x1)-func(t, T, e, x))/func(t, T, e, x))
        if relerr < error:
            return (x, func(t, T, e, x1), func(t, T, e, x),n)
        x = x1
        n += 1