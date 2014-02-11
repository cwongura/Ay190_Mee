# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 02:29:38 2014

@author: MeenoiKung
"""

import numpy as np
import time
import Ay190_Lib_9 as Lib9
import matplotlib.pyplot as pl
from plot_defaults import *

N = 5
a1 = np.loadtxt('LSE1_m.dat')
b1 = np.loadtxt('LSE1_bvec.dat')
a2 = np.loadtxt('LSE2_m.dat')
b2 = np.loadtxt('LSE2_bvec.dat')
a3 = np.loadtxt('LSE3_m.dat')
b3 = np.loadtxt('LSE3_bvec.dat')
a4 = np.loadtxt('LSE4_m.dat')
b4 = np.loadtxt('LSE4_bvec.dat')
a5 = np.loadtxt('LSE5_m.dat')
b5 = np.loadtxt('LSE5_bvec.dat')

t_Gauss = np.zeros(N)
t_Numpy = np.zeros(N)
a = np.array([a1, a2, a3, a4, a5])
b = np.array([b1, b2, b3, b4, b5])
n = np.array([len(b1),len(b2),len(b3),len(b4),len(b5)])

for i in range(N):
    t1 = time.clock()
    x = Lib9.GaussEli(a[i], b[i])
    t_Gauss[i] = time.clock()-t1

for i in range(N):
    t1 = time.clock()
    x = np.linalg.solve(a[i],b[i])
    t_Numpy[i] = time.clock()-t1

t_Gauss *= 1000
t_Numpy *= 1000

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

ax = pl.gca()
p1, = pl.plot(n, t_Gauss, 'b', linewidth=2)
p2, = pl.plot(n, t_Numpy, 'r', linewidth=2)

pl.legend((p1,p2),("Gaussian Elimination","LAPACK routine"),
          loc=(0.,0.7), frameon=False)
          
ax.set_xlabel('Size')
ax.set_ylabel('Time Taken')

ax.set_xscale('Log')
ax.set_yscale('Log')

pl.show()