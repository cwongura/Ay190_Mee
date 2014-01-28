# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 18:25:11 2014

Library for Ay190 Computational Astrophysics - WS 05

@author: Chatarin (Mee) Wong-u-railertkun
"""
#--------------------------------------------
import numpy as np
        
#--------------------------------------------
# Work Sheet 5
# Linear Regression
#--------------------------------------------

def LinRegress(x, y, sig):
    S = np.sum(1/sig**2)
    sumx = np.sum(x/sig**2)
    sumy = np.sum(y/sig**2)
    sumx2 = np.sum(x**2/sig**2)
    sumxy = np.sum(x*y/sig**2)
    
    a1 = (sumy*sumx2-sumx*sumxy)/(S*sumx2-sumx**2)
    a2 = (S*sumxy-sumy*sumx)/(S*sumx2-sumx**2)
    
    return (a1, a2)