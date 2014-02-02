# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 15:38:46 2014

@author: MeenoiKung
"""

def parabola(x):
    return x**2 + 1.

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_7 as Lib7

opt = 7
N_0 = 5
n = np.arange(opt)
n = N_0 * 10**(n)

ans = np.zeros(opt)
exact = 22/3.

for i in range(opt):
    ans[i] = Lib7.MCIntegrate(2, 3, parabola, n[i])

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(n, np.abs(exact-ans), "r", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of Points", labelpad = 15)
pl.ylabel("Alsolute Error", labelpad = -5)

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()