# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 07:17:32 2014

@author: MeenoiKung
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters
initial_data_file = "sun_earth.asc"
distance_unit_to_cm = 1.
time_unit_to_s = 1.
mass_unit_to_g = 1.

final_data_file = "final_positions_energy.asc"

def NbodyRHS(u,mass,time):
    x_old = u[:,0]
    y_old = u[:,1]
    z_old = u[:,2]
    
    n = len(x_old)    
    
    vx_old = u[:,3]
    vy_old = u[:,4]
    vz_old = u[:,5]
    vx_new = np.zeros(n)
    vy_new = np.zeros(n)
    vz_new = np.zeros(n)
    
    dx = vx_old
    dy = vy_old
    dz = vz_old
    
    for i in range(len(x_old)):
        xtemp = np.delete(x_old,i)
        ytemp = np.delete(y_old,i)
        ztemp = np.delete(z_old,i)
        masstemp = np.delete(mass,i)
        
        vec_x = x_old[i]-xtemp
        vec_y = y_old[i]-ytemp
        vec_z = z_old[i]-ztemp
        
        length = (vec_x**2+vec_y**2+vec_z**2)**(0.5)
        
        ax = -ggrav*masstemp*vec_x/length**2
        ay = -ggrav*masstemp*vec_y/length**2
        az = -ggrav*masstemp*vec_z/length**2

        vx_new[i] = np.sum(ax)
        vy_new[i] = np.sum(ay)
        vz_new[i] = np.sum(az)
    
    u_RHS = np.zeros(np.shape(u))
    u_RHS[:,0] = dx
    u_RHS[:,1] = dy
    u_RHS[:,2] = dz
    u_RHS[:,3] = vx_new
    u_RHS[:,4] = vy_new
    u_RHS[:,5] = vz_new
    return u_RHS

def NbodyRK4(u,mass,time,dt):
    k1 = dt*NbodyRHS(u, mass, time)
    k2 = dt*NbodyRHS(u+0.5*k1, mass, time+0.5*dt)
    k3 = dt*NbodyRHS(u+0.5*k2, mass, time+0.5*dt)
    k4 = dt*NbodyRHS(u+k3, mass, time+dt)
    
    y_new = u + (1./6)*(k1+2.*k2+2.*k3+k4)
    
    return y_new

def TotalEnergy(u,mass,time):
    x_old = u[:,0]
    y_old = u[:,1]
    z_old = u[:,2]
    vx_old = u[:,3]
    vy_old = u[:,4]
    vz_old = u[:,5]
    
    n = len(x_old)
    
    E_k = np.sum(0.5*mass*(vx_old**2+vy_old**2+vz_old**2))
    
    E_pot = 0.
    for i in range(n):
        xtemp = np.delete(x_old,i)
        ytemp = np.delete(y_old,i)
        ztemp = np.delete(z_old,i)
        masstemp = np.delete(mass,i)
        
        vec_x = x_old[i]-xtemp
        vec_y = y_old[i]-ytemp
        vec_z = z_old[i]-ztemp
        
        length = (vec_x**2+vec_y**2+vec_z**2)**(0.5)
        
        E_pot += np.sum(-ggrav*mass[i]*masstemp/length)
    
    return E_pot + E_k

# main program
plt.ion()

Nsteps = np.array([10,100,1000,10000])
ene_err = np.zeros(len(Nsteps))
for k in range(len(Nsteps)):
    t0 = 0
    t1 = 3 * seconds_per_year
    dt = (t1-t0)/Nsteps[k]

    (x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


    # convert from unitis in initial data file to cgs
    x *= distance_unit_to_cm
    y *= distance_unit_to_cm
    z *= distance_unit_to_cm
    vx *= distance_unit_to_cm / time_unit_to_s
    vy *= distance_unit_to_cm / time_unit_to_s
    vz *= distance_unit_to_cm / time_unit_to_s
    mass *= mass_unit_to_g

    xmin = np.amin(x)
    xmax = np.amax(x)
    ymin = np.amin(y)
    ymax = np.amax(y)
    zmin = np.amin(z)
    zmax = np.amax(z)
    rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

    # use a single state vector to simplify the ODE code
    # indices:
    # u[:,0] = x
    # u[:,1] = y
    # u[:,2] = z
    # u[:,3] = vx
    # u[:,4] = vy
    # u[:,5] = vz
    u = np.array((x,y,z,vx,vy,vz)).transpose()
    
    ini_ene = TotalEnergy(u, mass, 0)

    for it in range(0, Nsteps[k]):
        time = t0 + it * dt
        u = NbodyRK4(u,mass,time,dt)
        if it % max(1,Nsteps[k]/100) == 0:
            print "it = %d, time = %g years, energy = %g" % \
                (it, time / seconds_per_year,
                 TotalEnergy(u,mass,time))
            plt.clf()
            fig = plt.gcf()
            ax = mpl3d.Axes3D(fig)
            ax.scatter(u[:,0],u[:,1],u[:,2])
            ax.set_xlim((-rmax,rmax))
            ax.set_ylim((-rmax,rmax))
            ax.set_zlim((-rmax,rmax))
            plt.draw()
    
    ene_err[k] = np.abs(ini_ene-TotalEnergy(u,mass,time))

plt.figure()
plt.plot(Nsteps, ene_err)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Nsteps')
plt.ylabel('Change of energy')

plt.show()

# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)