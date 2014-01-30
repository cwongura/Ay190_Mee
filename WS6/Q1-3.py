# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:34:45 2014

Ay 190 - WS 06
Question 1b

@author: MeenoiKung
"""

import numpy as np
import Ay190_Lib_6 as Lib6
import matplotlib.pyplot as pl
from plot_defaults import *
from timeit import timeit

time1 = np.array([])
time2 = np.array([])
n = np.arange(10,1001,10)

for i in n:
    a =timeit("Lib6.dft(x)", number=10, \
              setup="import Ay190_Lib_6 as Lib6; import pylab; \
              x = pylab.randn(%d)" % i)
    b =timeit("np.fft.fft(x)", number=10, \
              setup="import numpy as np; import pylab; \
              x = pylab.randn(%d)" % i)
              
    time1 = np.append(time1,a)
    time2 = np.append(time2,b)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(n, time1, "r", linewidth=2)
p2, = pl.plot(n, time2, "b", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Length of input vector", labelpad = 15)
pl.ylabel("Computational time", labelpad = -5)

pl.legend((p1,p2),("dft(x)","fft(x)"),loc=(0.0,0.8),frameon=False)

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()