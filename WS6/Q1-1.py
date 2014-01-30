# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:44:41 2014

Ay 190 - WS 06
Question 1a

@author: MeenoiKung
"""

import numpy as np
import Ay190_Lib_6 as Lib6

n = 10
x = np.arange(n)

print "With input of (0,1,2,...,9)"
print "Result from dft(x) is:"
print Lib6.dft(x)
print "Result from fft(x) is:"
print np.fft.fft(x)