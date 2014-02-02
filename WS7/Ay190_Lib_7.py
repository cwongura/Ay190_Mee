# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:37:49 2014

Library for Ay190 Computational Astrophysics - WS 07

@author: Chatarin (Mee) Wong-u-railertkun
"""
#--------------------------------------------
import numpy as np

#--------------------------------------------
# Work Sheet 7
# Monte Carlo Experiment
#--------------------------------------------

def EstPi(n):
    np.random.seed(1)
    # Randomly choose the points, and move the poinst so the origin of the
    # circle is at the middle of the box
    x = np.random.rand(n) - 0.5
    y = np.random.rand(n) - 0.5
    
    m = 0 # Set the counter to zero
    r = x**2 + y**2
    
    for i in r:
        if i <= 0.5**2:
            m += 1
    return 4.*m/n

def EstPiSys(n):
    
    x = np.array([])
    y = np.array([])    
    
    import random
    sr = random.SystemRandom()
    for i in range(np.int(n)):
        x = np.append(x,sr.random() - 0.5)
        y = np.append(y,sr.random() - 0.5)
    
    m = 0
    r = x**2 + y**2
    
    for i in r:
        if i <= 0.5**2:
            m += 1
    return 4.*m/n

def Birthday(n):
    np.random.seed(1)
    
    trial = 1000
    m = 0
    for i in range(trial):
        ans = np.random.random_integers(1, 365, size=n)
        if len(np.unique(ans)) < n:
            m += 1
    return m*100/trial
    
def MCIntegrate(a, b, func, N):
    h = b - a
    fa = func(a)
    fb = func(b)
    
    A1 = fa*h
    A2 = (fb-fa)*h
    # Counter
    m = 0
    
    np.random.seed(1)
    x = np.random.uniform(a, b, size = N)
    y = np.random.uniform(fa, fb, size = N)
    for i in range(N):
        if func(x[i]) >= y[i]:
            m += 1
    return A2*m/N + A1
