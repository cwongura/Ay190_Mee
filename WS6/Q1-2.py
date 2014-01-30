# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:17:03 2014

Ay 190 - WS 06
Question 1b

@author: MeenoiKung
"""

import numpy as np
import Ay190_Lib_6 as Lib6
import matplotlib.pyplot as pl
from plot_defaults import *
from timeit import timeit

time = np.array([])
n = np.arange(10,1001,10)

for i in n:
    x =timeit("Lib6.dft(x)", number=10, \
              setup="import Ay190_Lib_6 as Lib6; import pylab; \
              x = pylab.randn(%d)" % i)
    time = np.append(time,x)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(n, time, "r", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Length of input vector", labelpad = 15)
pl.ylabel("Computational time", labelpad = -5)

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()