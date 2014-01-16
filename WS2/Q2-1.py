"""
Created on Mon Jan 13 21:29:08 2014

Ay 190 - WS 02
Question 2

@author: MeenoiKung
"""

import Ay190_Lib_2 as Lib2
import matplotlib.pyplot as pl
from plot_defaults import *

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

x1, y1 = Lib2.forward(Lib2.polyfunc, -3., 7., 1.)
x2, y2 = Lib2.central(Lib2.polyfunc, -3., 7., 1.)
y3 = Lib2.derifunc(x2)

p1, = pl.plot(x1, y1, "r", linewidth=2)
p2, = pl.plot(x2, y2, "b", linewidth=2)
p3, = pl.plot(x2, y3, "g", linewidth=2)

pl.axis([-2., 6., min(y1), max(y2)])

ax = pl.gca()

# label the axes
pl.xlabel("X",labelpad=15)
pl.ylabel("f'(X)",labelpad=-5)

pl.legend( (p1,p2,p3), ("Forward","Central","Analytic"), 
           loc=(0.29,0.5), frameon=False )
           
pl.savefig("CompareMethods.pdf")