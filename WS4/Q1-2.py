# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:56:15 2014

@author: MeenoiKung
"""

import numpy as np
import matplotlib.pyplot as pl
import Ay190_Lib_4 as Lib4

time = np.array([91., 182., 273.])
E_val = time * 0.
E_val = time * 0.
value_next = time * 0.
value = time * 0.
n = time * 0.

# Define the parameters we will use in this problem
period = 365.25635
e = 0.99999
a = 1.496 * 10**6
b = np.sqrt(a**2 * (1-e**2))

print "  t  ||        E      ||  Function Value    ||  x  ||  y  ||   n"
for i in range(3):
    (E_val[i], value_next[i], value[i], n[i]) = Lib4.Newton(Lib4.Efunc,
 Lib4.Edfunc, 0., 10**(-10), time[i], e = e)
    print time[i], "||", E_val[i], "||", value_next[i], "||", 
    print a*np.cos(E_val[i]), "||", b*np.sin(E_val[i]), "||", n[i], "\n"