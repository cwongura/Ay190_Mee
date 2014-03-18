# -*- coding: utf-8 -*-
"""
For Ay 190 - Computational Astrophysics - Final Project

Contain functions used to do numerical integration

@author: MeenoiKung
"""

import numpy as np

ggrav = 6.67e-8

def CalAcc(u, mass):
    '''
    Calculate acceleration of each particle in the system from gravitational
    force acted on the particle by the rest of the system.
    
    Let N be number of particle in the system
    Input: u - narray of shape (N, 3). Each column represents position of
               each particle in x, y, and z direction, respectively.
           mass - narray of length N. Contains mass of each particle.
    Output: narray of shape (N, 3). Each column represents accleration of 
            each particle in x, y, and z direction, respectively.
    '''
    x = u[:,0]
    y = u[:,1]
    z = u[:,2]
    
    n = len(x)
    
    a_x = np.zeros(n)
    a_y = np.zeros(n)
    a_z = np.zeros(n)
    
    for i in range(n):
        # Create temporary array that contains data for all particles in the
        # system, except the one in question.
        xtemp = np.delete(x,i)
        ytemp = np.delete(y,i)
        ztemp = np.delete(z,i)
        masstemp = np.delete(mass,i)
        
        # For each particle, find the vector pointing to all other particles
        # in the system
        vec_x = x[i]-xtemp
        vec_y = y[i]-ytemp
        vec_z = z[i]-ztemp
        
        # Calculate distance between the particle in question and all other
        # particles in the system
        length = (vec_x**2+vec_y**2+vec_z**2)**(0.5)
        
        # Calculate force exerted on the particle in question from all other
        # particles in the system
        fx = -ggrav*masstemp*vec_x/length**3
        fy = -ggrav*masstemp*vec_y/length**3
        fz = -ggrav*masstemp*vec_z/length**3
        
        # Sum up all the forces and divide by mass of the particle in question
        # in order to get the accleration
        a_x[i] = np.sum(fx)
        a_y[i] = np.sum(fy)
        a_z[i] = np.sum(fz)
    
    ans = np.zeros((n, 3))
    ans[:,0] = a_x
    ans[:,1] = a_y
    ans[:,2] = a_z
    return ans

def NbodyRHS(u,mass,time):
    '''
    Calculate the RHS of the system of equations,
    i.e. dx = velocidy and dv = acceleration
    
    Let N be number of particle in the system
    Input: u - narray of shape (N, 6). Each column represents position of
               each particle in x, y, and z direction, respectively, and 
               velocity of each particle in three directions.
           mass - narray of length N. Contains mass of each particle.
           time - float. The time of the system.
    Output: u_RHS - narray of shape (N, 6). First three columns are velocities
                    in x, y, z direction. The next three columns are
                    accelerations in x, y, z direction.
    '''
    vx_old = u[:,3]
    vy_old = u[:,4]
    vz_old = u[:,5]
    
    acceleration = CalAcc(u[:,0:3], mass)
    
    u_RHS = np.zeros(np.shape(u))
    u_RHS[:,0] = vx_old
    u_RHS[:,1] = vy_old
    u_RHS[:,2] = vz_old
    u_RHS[:,3:6] = acceleration

    return u_RHS
    
def NbodyRK4(u,mass,time, dt):
    '''
    Using RK4 numerical integration
    
    Let N be number of particle in the system
    Input: u - narray of shape (N, 6). Each column represents position of
               each particle in x, y, and z direction, respectively, and 
               velocity of each particle in three directions.
           mass - narray of length N. Contains mass of each particle.
           time - float. The time of the system.
           dt - size of time step
    Output: u_new - narray of shape (N, 6). Arranged in similar manner
                    as u, the array represents data at the next time step.
    '''
    k1 = dt*NbodyRHS(u, mass, time)
    k2 = dt*NbodyRHS(u+0.5*k1, mass, time+0.5*dt)
    k3 = dt*NbodyRHS(u+0.5*k2, mass, time+0.5*dt)
    k4 = dt*NbodyRHS(u+k3, mass, time+dt)
    
    y_new = u + (1./6)*(k1+2.*k2+2.*k3+k4)
    
    return y_new

def NbodyLeapFrog(u, mass, dt):
    '''
    Do the Leapfrog integration method which is
    x_{i+1} = x_i + v_i * dt + 0.5 * a_i * dt^2
    v_{i+1} = v_i + 0.5 * (a_i + a_{i+1}) * dt
    
    This is a sympletic method since the new velocity based on acceleration
    calculated from new position
    
    Let N be number of particle in the system
    Input: u - narray of shape (N, 6). Each column represents position of
               each particle in x, y, and z direction, respectively, and 
               velocity of each particle in three directions.
           mass - narray of length N. Contains mass of each particle.
           dt - size of time step
    Output: u_new - narray of shape (N, 6). Arranged in similar manner
                    as u, the array represents data at the next time step.
    '''
    x_old = u[:,0]
    y_old = u[:,1]
    z_old = u[:,2]
    vx_old = u[:,3]
    vy_old = u[:,4]
    vz_old = u[:,5]
    
    u_new = np.zeros(np.shape(u))
    
    # Calculate acceleration based on the position now
    acc_now = CalAcc(u[:,0:3],mass)
    acc_nx = acc_now[:,0]
    acc_ny = acc_now[:,1]
    acc_nz = acc_now[:,2]
    
    # Calculate the new position 
    x_new = x_old + vx_old*dt + 0.5*acc_nx*dt*dt
    y_new = y_old + vy_old*dt + 0.5*acc_ny*dt*dt
    z_new = z_old + vz_old*dt + 0.5*acc_nz*dt*dt
    
    # Assign the new position into the new u array
    u_new[:,0] = x_new
    u_new[:,1] = y_new
    u_new[:,2] = z_new
    
    # Calculate acceleration based on the new position
    acc_future = CalAcc(u_new[:,0:3],mass)
    acc_fx = acc_future[:,0]
    acc_fy = acc_future[:,1]
    acc_fz = acc_future[:,2]
    
    # Calculate the new velocity
    vx_new = vx_old + 0.5*(acc_nx+acc_fx)*dt
    vy_new = vy_old + 0.5*(acc_ny+acc_fy)*dt
    vz_new = vz_old + 0.5*(acc_nz+acc_fz)*dt
    
    # Assign the new velocity into the new u array
    u_new[:,3] = vx_new
    u_new[:,4] = vy_new
    u_new[:,5] = vz_new
    
    return u_new