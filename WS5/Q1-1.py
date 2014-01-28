"""
Created on Tue Jan 21 09:50:33 2014

For Ay190 - WS 05

@author: MeenoiKung
"""

#!/usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *

data = ascii.read("m_sigma_table.dat",readme="m_sigma_ReadMe.dat")
logsigma = np.array(np.log10(data["sigma*"]))
logM = np.array(data["logM"])
d_logM = np.array(data["e_logM"])

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

(a, b, c) = pl.errorbar(logsigma, logM, yerr=d_logM, fmt='o')

ax = pl.gca()

# Label the axes
pl.xlabel("log($\sigma_*$ / km $s^{-1}$)",labelpad=15)
pl.ylabel("log($M_{BH}/M_{solar}$)",labelpad=5)

pl.show()