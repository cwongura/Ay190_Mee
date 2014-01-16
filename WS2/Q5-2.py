# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 07:05:14 2014

Ay 190 - WS 02
Question 4

@author: MeenoiKung
"""

import Ay190_Lib_2 as Lib2
import matplotlib.pyplot as pl
import scipy.interpolate
from plot_defaults import *

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

node = np.array([0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
val = np.array([0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561,
                    0.468, 0.302])

x1, y1 = Lib2.HermiteInter(node,val)

tck = scipy.interpolate.splrep(node,val,s=0)
xnew = np.linspace(0.,1.,num=10**3)
ynew = scipy.interpolate.splev(xnew,tck,der=0)

p1, = pl.plot(node,val,'ro',markersize=17.)
p2, = pl.plot(x1, y1, 'b', linewidth=2)
p3, = pl.plot(xnew, ynew, 'g', linewidth=2)

ax = pl.gca()

# label the axes
pl.xlabel("Time (days)",labelpad=15)
pl.ylabel("Apparent Magnitude",labelpad=5)

pl.legend( (p2,p3), ("$Hermite$","$Spline$"), 
           loc=(0.6,0.15), frameon=False )

pl.show()
pl.savefig("SplineInter.pdf")