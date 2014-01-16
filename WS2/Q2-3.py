# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 01:02:48 2014

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

x1, y1 = Lib2.central(Lib2.polyfunc, -2., 6., 10**(-3))
x2, y2 = Lib2.central(Lib2.polyfunc, -2., 6., (0.5*10**(-3)))

error1 = Lib2.derifunc(x1)-y1
error2 = Lib2.derifunc(x2)-y2

p1, = pl.plot(x1,error1, "r", linewidth=2)
p2, = pl.plot(x2,error2, "b", linewidth=2)


ax = pl.gca()

# label the axes
pl.xlabel("X",labelpad=15)
pl.ylabel("Absolute Error",labelpad=-8)

pl.legend( (p1,p2), ("$h_1$","$h_2=h_1/2$"), 
           loc=(0.4,0.3), frameon=False )

pl.show()
pl.savefig("CentralErr.pdf")