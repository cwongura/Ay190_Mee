# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 09:36:32 2014

@author: MeenoiKung
"""

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
import Ay190_Lib_7 as Lib7

n = np.arange(2, 50)
prob = np.zeros(len(n))
half = prob+50
for i in n:
    prob[i-2] = Lib7.Birthday(i)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

p1, = pl.plot(n, prob, "r", linewidth=2)
p2, = pl.plot(n, half, "b", linewidth=2)

ax = pl.gca()

# Label the axes
pl.xlabel("Number of Student", labelpad = 15)
pl.ylabel("Probability", labelpad = -5)

print prob

pl.show()