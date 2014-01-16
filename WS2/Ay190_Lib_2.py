"""
Created on Mon Jan 13 21:29:08 2014

Library for Ay190 Computational Astrophysics - WS 02

@author: MeenoiKung
"""
import numpy as np

#--------------------------------------------
# Work Sheet 2
#--------------------------------------------

#--------------------------------------------
# Question 1
#--------------------------------------------
def recur13(n):
    '''
    Return the nth order term of (1/3)^n using single precision
    FP numbers (float32)
    '''
    n = int(n) # Make sure that n is an integer
    if n < 0: # If n is negative, tell the user and return 0.
        print "Error: n has to be positive"
        return 0
        
    ans = np.array([1., 1./3], dtype = np.float32)
    if n == 0 or n == 1:
        return ans[n]
    else:
        for i in range(2, n+1):
            ans = np.append(ans, 0.)
            ans[i] = np.float32((13.*ans[i-1]/3.) - (4.*ans[i-2]/3.))
        return np.float32(ans[n])

#--------------------------------------------
# Question 2
#--------------------------------------------

def polyfunc(x):
    '''
    Return the function x^3 - 5x^2 + x
    '''
    return x**3 - 5.* (x**2) + x

def derifunc(x):
    '''
    Return the function 3*x^2 - 10*x + 1
    '''
    
    return 3. * x**2 - 10. * x + 1.
    
def central(func, a, b, h):
    node = np.linspace(a, b, num = (b-a)/h)
    ans = node*0.
    
    for i in range(1, len(node)-1):
        ans[i] = (func(node[i+1]) - func(node[i-1]))/(2.*h)
    
    return node[1:-1], ans[1:-1]
    
def forward(func, a, b, h):
    node = np.linspace(a, b, num = (b-a)/h)
    ans = node*0.
    
    for i in range(len(node)-1):
        ans[i] = (func(node[i+1]) - func(node[i]))/h
    return node[:-1], ans[:-1]
    
#--------------------------------------------
# Question 4
#--------------------------------------------
def LagrangeInter(x, node, val):                
    ans = 0.
    for i in range(len(node)):
        x_j = node[i]
        m = 1.
        for x_k in node:
            if x_j == x_k:
                pass
            else:
                m *= (x-x_k)/(x_j - x_k)
        ans += val[i] * m
    return ans

def LinearInter(node, val):
    x = np.array([])
    ans = np.array([])
    
    for i in range(len(node)-1):
        space = np.linspace(node[i],node[i+1],num = (node[i+1]-node[i])*10**2)
        x = np.append(x, space)
        ans = np.append(ans, val[i]+(val[i+1]-val[i])/(node[i+1]-node[i])
        *(space-node[i]))
    return x, ans

def QuadInter(node, val):
    x = np.array([])
    ans = np.array([])
    
    for i in range(len(node)-2):
        space = np.linspace(node[i],node[i+1],num = (node[i+1]-node[i])*10**2)
        x = np.append(x, space)
        ans = np.append(ans,(space-node[i+1])*(space-node[i+2])*val[i]/
        (node[i]-node[i+1])/(node[i]-node[i+2])+(space-node[i])*(space-node[i+2])
        *val[i+1]/(node[i+1]-node[i])/(node[i+1]-node[i+2])+(space-node[i])*
        (space-node[i+1])*val[i+2]/(node[i+2]-node[i])/(node[i+2]-node[i+1]))
    return x, ans

#--------------------------------------------
# Question 5
#--------------------------------------------
def psi0(z):
    return 2.*z**3 - 3.*z**2 + 1
    
def psi1(z):
    return z**3 - 2.*z**2 + z
    
def HermiteInter(node, val):
    deri = (val[1:] - val[:-1])/(node[1:] - node[:-1])
    x = np.array([])
    ans = np.array([])
    
    for i in range(len(node)-2):
        space = np.linspace(node[i],node[i+1], num = (node[i+1]-node[i])*10**2)
        x = np.append(x, space)
        z = (space-node[i])/(node[i+1]-node[i])
        
        ans = np.append(ans,val[i]*psi0(z)+val[i+1]*psi0(1-z)+deri[i]*
        (node[i+1]-node[i])*psi1(z)-deri[i+1]*(node[i+1]-node[i])*psi1(1-z))
    return x, ans