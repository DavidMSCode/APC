# %%
#External functions
#import modules
import sys, os
import APC as APC
import numpy as np
from Kepler import elms2rv, muEarth
import matplotlib.pyplot as plt
import pandas as pd
# import faulthandler
print('finished')

# %%
q = 8000
e = .05
a = q/(1-e)
i = 63.435*np.pi/180
r0,v0 = elms2rv(a,e,i,0,0,0,muEarth)
T = 2*np.pi*np.sqrt(a**3/muEarth)
t0 = 0.0
tf = 1000*T
#LEO
# r0 = [6500, 0.0, 0.0];                                    # Initial Position (km)
# v0 = (0, 7.8309, 0.0);                              # Initial Velocity (km/s)
# t0 = 0.0;                                                 # Initial Times (s)
# T = np.pi*np.sqrt(6500*3/muEarth)
# tf = 100*T;                                     # Final Time (s)
# #MEO
# r0 = [9000.0, 0.0, 0.0];                                    # Initial Position (km)
# v0 = [0.0, 6.7419845635570, 1.806509319188210];             # Initial Velocity (km/s)
# t0 = 0.0;                                                   # Initial Times (s)
# tf = 10*9.952014050491189e+03;                              # Final Time (s)

#sat props
mass = 1000;
area = 10;
reflectance = 1.5;
Cd = 2.0;

# #run APC code in single orbit mode with and without perturbations
#output = APC.SinglePropagate(r0,v0,t0,tf,area,reflectance,mass,Cd,False,False,False)
# #use a user specified time vec
time_vec = np.arange(t0,tf+1,300)
output2 = APC.SinglePropagate(r0,v0,time_vec,area,reflectance,mass,Cd,False,False,False)
ax = plt.axes(projection='3d')
ax.plot3D(output2.getPositionX(),output2.getPositionY(),output2.getPositionZ())
plt.show()

# %%
