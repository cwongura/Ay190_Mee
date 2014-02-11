# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 01:07:49 2014

For Ay 190 - Computational Astrophysics - Library WS 09

@author: MeenoiKung
"""

import numpy as np
from scipy.linalg import lu

def GaussEli(a, b):
    N = len(b)
    x = np.zeros(N)
    
    (pl, u) = lu(a, permute_l=True)
    pl = np.mat(pl)

    bvec = np.reshape(np.mat(b),(N,1))
    b = np.linalg.inv(pl)*bvec
   
    for i in range(N-1, -1, -1):
        ans = b[i]/u[i,i]

        for n in range(i, N-1):
            ans -= u[i,n+1]*x[n+1]/u[i,i]
        x[i] = ans
        
    return x