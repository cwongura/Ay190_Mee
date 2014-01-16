"""
Created on Thu Jan 16 02:01:37 2014

Ay 190 - WS 02
Question 4

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

node = np.array([0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
val = np.array([0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561,
                    0.468, 0.302])
x = np.linspace(0., 1., num = 10**3)

p1, = pl.plot(node,val,'ro',markersize=17.)
p2, = pl.plot(x,Lib2.LagrangeInter(x, node, val), 'b', linewidth=2)

ax = pl.gca()

# label the axes
pl.xlabel("Time (days)",labelpad=15)
pl.ylabel("Apparent Magnitude",labelpad=5)

#pl.legend( (p1,p2), ("$h_1$","$h_2=h_1/2$"), 
#           loc=(0.6,0.75), frameon=False )

pl.show()
pl.savefig("GlobalLagrange.pdf")