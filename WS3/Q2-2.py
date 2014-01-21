# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 01:04:58 2014

For Ay190 - WS 04

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_3 as Lib3
from plot_defaults import *

beta = 1/20.
h_bar = 6.582119*10**(-22)
c = 3*10**8

bins = np.arange(0, 151, 5)
N = np.array([1,10,100,10**3,10**4])
ans = N * 0.

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

for i in range(len(N)-1):
    total = 0.
    for m in range(len(bins)-1):
        total += Lib3.GaussLegendre(Lib3.EnergyFunc, bins[m], bins[m+1],N[i])
    ans[i] = total

ans *= (8*np.pi)*5./(2*np.pi*h_bar*c)**3/10**42

p1, = pl.plot(N, ans, "r", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of nodes", labelpad = 15)
pl.ylabel("Total number density ($10^{42} m^{-3}$)", labelpad = 0)

pl.show()
pl.savefig("EnergyDens.pdf")