# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 01:27:08 2014

@author: MeenoiKung
"""

import numpy as np
import Ay190_Lib_9 as Lib9

np.random.seed(8)
a = np.reshape(np.random.rand(9),(3,3))
x = np.reshape(np.random.rand(3),(3,1))
b = np.reshape(np.asarray(np.mat(a)*np.mat(x)),3)
x = np.reshape(x,3)

x_Gauss = Lib9.GaussEli(a, b)

print "The absolute error is", np.abs(x-x_Gauss)