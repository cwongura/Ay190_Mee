# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 04:24:02 2014

For Ay190 - WS 03

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_3 as Lib3
from plot_defaults import *

a = 0.
b = np.pi
N = np.arange(10, 10**3)
mid = N * 0.
trap = N * 0.
simp = N * 0.
mid2 = N * 0.

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

true_val = 2.

for i in range(len(N)):
    mid[i] = abs(true_val-Lib3.Midpoint(np.sin, a, b, N[i]))
    simp[i] = abs(true_val-Lib3.Simpson(np.sin, a, b, N[i]))
    trap[i] = abs(true_val-Lib3.Trapezoid(np.sin, a, b, N[i]))

p1, = pl.plot(N, mid, "r", linewidth=2)
p2, = pl.plot(N, trap, "b", linewidth=2)
p3, = pl.plot(N, simp, "g", linewidth=2)

ax = pl.gca()

# label the axes
pl.xlabel("Number of subsection",labelpad=15)
pl.ylabel("Absolute Error",labelpad=-8)

ax.set_xscale('log')
ax.set_yscale('log')

pl.legend((p1, p2, p3), ("Midpoint", "Trapezoid", "Simpson"),
          loc=(0.6, 0.3), frameon = False)

pl.show()
pl.savefig("ConvergentA.pdf")