"""
Created on Mon Jan 13 21:29:08 2014

Library for Ay190 Computational Astrophysics - WS 03

@author: Chatarin (Mee) Wong-u-railertkun
"""
#--------------------------------------------
import numpy as np
from scipy import special as sp
        
#--------------------------------------------
# Work Sheet 3
# Numerical Integration
#--------------------------------------------

#--------------------------------------------
# Question 1
#--------------------------------------------

def XSin(x):
    return x*np.sin(x)
        
def Midpoint(func, a, b, N):
    '''
    Numerical integration using midpoint method
    i.e., Q_i = h * f((a_i + b_i)/2)
    
    Input: func - the integrand function
           a    - initial point
           b    - final point
           N    - number of quadrature
    '''
    
    node = np.linspace(a, b, num = N+1)
    mid = 0.5 * (node[1:] + node[:-1])
    h = node[1]-node[0]
    
    return np.sum(h * (func(mid)))

def Trapezoid(func, a, b, N):
    '''
    Numerical integration using the trapezoid rule
    i.e., Q_i = 0.5 * h * (f(b_i) + f(a_i))
    
    Input: func - the integrand function
           a    - initial point
           b    - final point
           N    - number of quadrature
    '''
    node = np.linspace(a, b, num = N+1)
    h = node[1]-node[0]
    
    return h * (np.sum(func(node[1:-1])) + 0.5*(func(node[0])+func(node[-1])))
    
def Simpson(func, a, b, N):
    '''
    Numerical integration using the Simpson's rule
    i.e., Q_i = (h/6)*(f(a_i)+4*f((a_i+b_i)/2)+f(b_i))
    
    Input: func - the integrand function
           a    - initial point
           b    - final point
           N    - number of spacing
    '''
    if N%2==1:
        N += 1
    
    node = np.linspace(a, b, num = N+1)
    h = node[1]-node[0]
    mid = 0.5*(node[1:]+node[:-1])
    
    return h*np.sum(func(node[1:])+func(node[:-1])+4*func(mid))/6.

#--------------------------------------------
# Question 2
#--------------------------------------------

def ExFunc(x):
    return np.exp(x)*x**2/(np.exp(x)+1)

def EnergyFunc(E):
    beta = 1/20.
    return E**2/(np.exp(beta*E)+1)

def GaussLaguerre(func, n):
    [laguerre_roots, laguerre_weights] = sp.l_roots(n, 0.)
    return np.sum(laguerre_weights*func(laguerre_roots))

def GaussLegendre(func, a, b, n):
    [legendre_roots, legendre_weights] = sp.p_roots(n, 0.)
    
    return np.sum(legendre_weights*(b-a)*func((b-a)*
    legendre_roots/2.+(a+b)/2.)/2.)