# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 00:10:54 2014

For Ay 190 - Computational Astrophysics - WS 09

@author: MeenoiKung
"""

import numpy as np

a1 = np.loadtxt('LSE1_m.dat')
b1 = np.loadtxt('LSE1_bvec.dat')
a2 = np.loadtxt('LSE2_m.dat')
b2 = np.loadtxt('LSE2_bvec.dat')
a3 = np.loadtxt('LSE3_m.dat')
b3 = np.loadtxt('LSE3_bvec.dat')
a4 = np.loadtxt('LSE4_m.dat')
b4 = np.loadtxt('LSE4_bvec.dat')
a5 = np.loadtxt('LSE5_m.dat')
b5 = np.loadtxt('LSE5_bvec.dat')

(sign1, logdet1) = np.linalg.slogdet(a1)
(sign2, logdet2) = np.linalg.slogdet(a2)
(sign3, logdet3) = np.linalg.slogdet(a3)
(sign4, logdet4) = np.linalg.slogdet(a4)
(sign5, logdet5) = np.linalg.slogdet(a5)

print "  i  ,  shape  ,  determinant  ,  b size  "
print "1", np.shape(a1), np.linalg.det(a1), np.shape(b1)
print "2", np.shape(a2), np.linalg.det(a2), np.shape(b2)
print "3", np.shape(a3), np.linalg.det(a3), np.shape(b3)
print "4", np.shape(a4), np.linalg.det(a4), np.shape(b4)
print "5", np.shape(a5), np.linalg.det(a5), np.shape(b5)


