# %%
#External functions
#import modules
import sys, os
from time import time
from turtle import xcor
import APC as APC
import numpy as np
from Kepler import elms2rv, muEarth, muMoon
import matplotlib.pyplot as plt
import pandas as pd
# import faulthandler
print('finished')

# %%
Re = 1.7380000000000000e+03
q = Re+100
e = .05
a = q/(1-e)
i = 63.435*np.pi/180
r0,v0 = elms2rv(a,e,i,0,0,0,muMoon)
T = 2*np.pi*np.sqrt(a**3/muMoon)
t0 = 0.0
tf = 2*T

#sat props
mass = 1000;
area = 10;
reflectance = 1.5;
Cd = 2.0;
# %%
# #run APC code in single orbit mode with and without perturbations
output = APC.SinglePropagate(r0,v0,t0,tf,area,reflectance,mass,Cd,compute_hamiltonian=True,compute_drag=False)
# #use a user specified time vec
tstep = T/12
#time vec with higher density as t approaches tf
time_vec = np.flip(tf-np.geomspace(1,tf,32))
output2 = APC.SinglePropagate(r0,v0,time_vec,area,reflectance,mass,Cd,compute_hamiltonian=True)


#Plotting
l1 = "default 30s timestep"
l2 = "User specified geometric series timestep".format(tstep)
fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(1,2,1,projection='3d')
ax1.plot3D(output.getPositionX(),output.getPositionY(),output.getPositionZ(),label=l1)
ax1.plot3D(output2.getPositionX(),output2.getPositionY(),output2.getPositionZ(),'s',label=l2)
ax1.set_xlabel('X (km))')
ax1.set_ylabel('Y (km)')
ax1.set_zlabel('Z (km)')
ax2 = fig.add_subplot(1,2,2)
ax2.plot(output.getTimes()/T,output.getHamiltonian(),'-',label=l1)
ax2.plot(output2.getTimes()/T,output2.getHamiltonian(),'s',markersize=6,label=l2)
ax2.grid()
ax2.set_xlabel("Orbit #",fontsize=14)
ax2.set_ylabel(r"Hamiltonian $\frac{E-E_0}{E_0}$",fontsize=14)
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=1, 
                    top=1, 
                    wspace=0.6, 
                    hspace=0.4)
plt.legend()
plt.show()

# %%
# #use a user specified time vec
tstep = T/12
time_vec = np.arange(t0,tf+1,tstep)
output2 = APC.SinglePropagate(r0,v0,time_vec,area,reflectance,mass,Cd,compute_hamiltonian=True)


#Plotting
l1 = "default 30s timestep"
l2 = "User specified {:.0f}s timestep".format(tstep)
fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(1,2,1,projection='3d')
ax1.plot3D(output.getPositionX(),output.getPositionY(),output.getPositionZ(),label=l1)
ax1.plot3D(output2.getPositionX(),output2.getPositionY(),output2.getPositionZ(),'s',label=l2)
ax1.set_xlabel('X (km))')
ax1.set_ylabel('Y (km)')
ax1.set_zlabel('Z (km)')
ax2 = fig.add_subplot(1,2,2)
ax2.plot(output.getTimes()/T,output.getHamiltonian(),'-',label=l1)
ax2.plot(output2.getTimes()/T,output2.getHamiltonian(),'s',markersize=6,label=l2)
ax2.grid()
ax2.set_xlabel("Orbit #",fontsize=14)
ax2.set_ylabel(r"Hamiltonian $\frac{E-E_0}{E_0}$",fontsize=14)
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=1, 
                    top=1, 
                    wspace=0.6, 
                    hspace=0.4)
plt.legend()
plt.show()
# %%
