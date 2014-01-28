# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:51:29 2014

For Ay190 - WS 05

@author: MeenoiKung
"""

#!/usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_5 as Lib5

data = ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat")
logsigma = np.array(np.log10(data["sigma*"]))
logM = np.array(data["logM"])
d_logsigma = np.array(data["e_sigma*"])
d_logM = np.array(data["e_logM"])

(a1, a2) = Lib5.LinRegress(logsigma, logM, np.ones(len(logM)))
# We will use a2 as an approximate slope dy/dx

sig = np.sqrt(d_logM**2+(a2*d_logsigma)**2)

(a3, a4) = Lib5.LinRegress(logsigma, logM, sig)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(logsigma, logM, "ro", markersize=10)
p2, = pl.plot(logsigma, a1+a2*logsigma, "b", linewidth=2)
p3, = pl.plot(logsigma, a3+a4*logsigma, "g", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("log($\sigma_*$ / km $s^{-1}$)",labelpad=15)
pl.ylabel("log($M_{BH}/M_{solar}$)",labelpad=5)

# Legend
pl.legend((p2, p3), ("Without Errors","With Errors"),
          loc=(0.5,0.1),frameon=False)

print a3
print a4

pl.show()