# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:44:56 2014

For Ay190 - WS 07

@author: MeenoiKung
"""

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_7 as Lib7

m = 5
n = np.arange(0,m,0.25)
n = 10**n
#n = np.arange(1,10**3,10)

PiRand = np.zeros(len(n))
PiSysRand = np.zeros(len(n))

for i in range(len(n)):
    PiRand[i] = np.abs(np.pi - Lib7.EstPi(n[i]))
    PiSysRand[i] = np.abs(np.pi - Lib7.EstPiSys(n[i]))

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(n, PiRand, "r", linewidth=2)
p2, = pl.plot(n, PiSysRand, "b", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of random points", labelpad = 15)
pl.ylabel("Absolute Error", labelpad = -5)

# Legend
pl.legend((p1, p2), ("NumPy Random","System Random"),
          loc = (0.0,0.0), frameon=False)

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()