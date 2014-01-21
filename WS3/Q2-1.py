# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 06:47:53 2014

For Ay190 - WS 04

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_3 as Lib3
from plot_defaults import *

cons = 8.4383
N = np.arange(1, 15)
ans = N * 0.

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

for i in range(len(N)):
    ans[i] = cons*Lib3.GaussLaguerre(Lib3.ExFunc, N[i])

p1, = pl.plot(N, ans, "r", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of nodes", labelpad = 15)
pl.ylabel("Total number density ($10^{33} m^{-3}$)", labelpad = 0)

#ax.set_xscale('log')
#ax.set_yscale('log')

pl.show()
pl.savefig("NumDens.pdf")