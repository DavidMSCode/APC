# %%
#External functions
from draw_funcs import *
from Kepler import *
from Utils import out2breaks
#import modules
import sys, os, platform
import APC as APC
import plotly.graph_objects as go
import numpy as np
from IPython.display import display, Math, Latex
import spiceypy.spiceypy as sp
import matplotlib.pyplot as plt
import random

#%%
sp.furnsh("de421.bsp")
sp.furnsh("grail_120301_120529_sci_v02.bsp")
sp.furnsh("naif0012.tls")
mass = 212;
area = 10;
reflectance = 1.5;
Cd = 2.0;

# %%
int = 2*60*60
randoffset = 2*int
propagation_length = 24*60*60
date1 = "2012 April 1 12:00:00"
date2 = "2012 April 15 12:00:00"
et0 = sp.str2et(date1)
etf = sp.str2et(date2)
ts1 = np.arange(et0,etf,int)
ts1 = ts1 - np.random.rand(len(ts1))*randoffset-randoffset/2   #add random noise to start of propagation so that they aren't all near each other

# %%
errs = []
for t0 in ts1:
    tf = t0+propagation_length
    state0 = sp.spkezr("GRAIL-A",t0,"J2000","LT+S","MOON")[0]
    title = "Propagated lunar orbit"
    r0 = state0[0:3]
    v0 = state0[3:]
    #propagate initial state
    orbit1 = APC.Orbit("MOON","MOON_PA","J2000")
    orbit1.SetPosition0(r0)
    orbit1.SetVelocity0(v0)
    orbit1.SetIntegrationTime(t0,tf)
    orbit1.SetComputeThirdBody(True)
    orbit1.SetComputeSRP(False)
    orbit1.SetComputeHamiltonian(False)
    orbit1.SinglePropagate()
    print(".")
    #Get ephemeris position
    ts = np.array(orbit1.getTimes())+t0
    prop_pos = np.array(orbit1.getPosition()).T
    ephem_pos = sp.spkpos("GRAIL-A",ts,"J2000","LT+S","MOON")[0]
    err = prop_pos-ephem_pos
    errs.append(err)

# %%
path="images/"
title = "GRAIL-A propagation drift over {:.1f} hours for random initial states between ".format((propagation_length/3600)) +date1 + " & " + date2
ts = np.array(orbit1.getTimes())/3600
fig1 = plt.figure(figsize=(12,6))
ax1 = fig1.gca()
fig2 = plt.figure(figsize=(12,6))
ax2 = fig2.gca()
fig3 = plt.figure(figsize=(12,6))
ax3 = fig3.gca()
outlier=200
for i in reversed(range(0,len(errs))):
    err = errs[i]
    if (abs(errs[i])>outlier).any():
        #remove from list if outlier
        del errs[i]
        print("found outlier at %i",i)
    else:
        #otherwise plot
        ax1.plot(ts,err[:,0],color='black')
        ax2.plot(ts,err[:,1],color='black')
        ax3.plot(ts,err[:,2],color='black')

aerrs = np.sum(np.array(errs),axis=0)/len(errs)
ax1.plot(ts,aerrs[:,0],color='red',label="Avg")
ax2.plot(ts,aerrs[:,1],color='red',label="Avg")
ax3.plot(ts,aerrs[:,2],color='red',label="Avg")
ax1.grid()
ax2.grid()
ax3.grid()
ax1.set_xlabel("Time (hr)")
ax1.set_ylabel("J2000 X Error (km)")
ax1.set_title(title)
ax1.legend()
ax2.set_xlabel("Time (hr)")
ax2.set_ylabel("J2000 Y Error (km)")
ax2.set_title(title)
ax2.legend()
ax3.set_xlabel("Time (hr)")
ax3.set_ylabel("J2000 Z Error (km)")
ax3.set_title(title)
ax3.legend()
fig1.savefig(path+"X-Drift "+date2+".png",dpi=300,bbox_inches="tight")
fig2.savefig(path+"Y-Drift "+date2+".png",dpi=300,bbox_inches="tight")
fig3.savefig(path+"Z-Drift "+date2+".png",dpi=300,bbox_inches="tight")

# %%
# #run APC code in single orbit mode with and without perturbations
# get the best orbit
title="Example error for 24 hour GRAIL-A Propagation"
best = 0
oldavg = 10000000000
for i in range(1,len(errs)):
    avgerror = np.mean(np.linalg.norm(errs[i],axis=1))
    if avgerror<oldavg:
        oldavg=avgerror
        best = i
print("the best orbit is #{:d}".format(best))
fig1 = plt.figure(figsize=(12,6))
ax1 = fig1.gca()
fig2 = plt.figure(figsize=(12,6))
ax2 = fig2.gca()
fig3 = plt.figure(figsize=(12,6))
ax3 = fig3.gca()
fig4 = plt.figure(figsize=(12,6))
ax4 = fig4.gca()
err = errs[best]
ax1.plot(ts,err[:,0],color='black')
ax2.plot(ts,err[:,1],color='black')
ax3.plot(ts,err[:,2],color='black')
ax4.plot(ts,np.linalg.norm(err,axis=1))
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax1.set_xlabel("Time (hr)")
ax1.set_ylabel("J2000 X Error (km)")
ax1.set_title(title)
ax1.legend()
ax2.set_xlabel("Time (hr)")
ax2.set_ylabel("J2000 Y Error (km)")
ax2.set_title(title)
ax2.legend()
ax3.set_xlabel("Time (hr)")
ax3.set_ylabel("J2000 Z Error (km)")
ax3.set_title(title)
ax3.legend()
ax4.set_xlabel("Time (hr)")
ax4.set_ylabel("J2000 Total Error (km)")
ax4.set_title(title)
ax4.legend()
fig1.savefig(path+"X-Drift best.png",dpi=300,bbox_inches="tight")
fig2.savefig(path+"Y-Drift best.png",dpi=300,bbox_inches="tight")
fig3.savefig(path+"Z-Drift best.png",dpi=300,bbox_inches="tight")
fig4.savefig(path+"Total-Drift best.png",dpi=300,bbox_inches="tight")
# %%
ts = np.array(orbit1.getTimes())+et0
prop_pos = np.array(orbit1.getPosition()).T
ephem_pos = sp.spkpos("GRAIL-A",ts,"J2000","LT+S","MOON")[0]
np.linalg.norm(prop_pos-ephem_pos,axis=1)
# %%
date1 = "2012 April 1 12:00:00"
date2 = "2012 April 15 12:00:00"
et0 = sp.str2et(date1)
etf = sp.str2et(date2)
ts1 = np.arange(et0,etf,int)
ts1 = ts1 - np.random.rand(len(ts1))*randoffset-randoffset/2   #add random noise to start of propagation so that they aren't all near each other

# %%
tf = t0+5*24*3600
state0 = sp.spkezr("GRAIL-A",t0,"J2000","LT+S","MOON")[0]
title = "Propagated lunar orbit"
r0 = state0[0:3]
v0 = state0[3:]
#propagate initial state
orbit1 = APC.Orbit("MOON","MOON_PA","J2000")
orbit1.SetPosition0(r0)
orbit1.SetVelocity0(v0)
orbit1.SetIntegrationTime(t0,tf)
orbit1.SetComputeThirdBody(True)
orbit1.SetComputeSRP(False)
orbit1.SetComputeHamiltonian(False)
orbit1.SinglePropagate()
# %%
ts = np.array(orbit1.getTimes())+t0
ephem_pos = sp.spkpos("GRAIL-A",ts,"J2000","LT+S","MOON")[0]
#plot the full orbit
fig = plotOrbit(orbit1)
fig.update_layout(title=title)
fig.show()

#plot a closeup of the intial position
camera = dict(
    up=dict(x=0, y=0, z=1),
    eye=dict(x=.08, y=-.1, z=.05),
    center=dict(x=0,y=0,z=0)
)
diff = 6378
x0 = r0[0]
y0 = r0[1]
z0 = r0[2]
fig.update_layout(
    scene=dict(
    aspectmode = "cube",
    xaxis=dict(range=[x0-diff,x0+diff]),
    yaxis=dict(range=[y0-diff,y0+diff]),
    zaxis=dict(range=[z0-diff,z0+diff])),
    scene_camera=camera, 
    title="Propagated Position Close-up")
fig.show()

fig = plotEphem(ephem_pos)
fig.update_layout(title="Ephemeris Orbit")
fig.show()

#plot a closeup of the intial position
camera = dict(
    up=dict(x=0, y=0, z=1),
    eye=dict(x=.08, y=-.1, z=.05),
    center=dict(x=0,y=0,z=0)
)
diff = 6378
x0 = r0[0]
y0 = r0[1]
z0 = r0[2]
fig.update_layout(
    scene=dict(
    aspectmode = "cube",
    xaxis=dict(range=[x0-diff,x0+diff]),
    yaxis=dict(range=[y0-diff,y0+diff]),
    zaxis=dict(range=[z0-diff,z0+diff])),
    scene_camera=camera, 
    title="Ephemeris Position Close-up")
fig.show()
# %%
t0 = et0
tf = t0+propagation_length
state0 = sp.spkezr("GRAIL-A",t0,"J2000","LT+S","MOON")[0]
title = "Propagated lunar orbit"
r0 = state0[0:3]
v0 = state0[3:]
#propagate initial state
orbit1 = APC.Orbit("MOON","MOON_PA","J2000")
orbit1.SetPosition0(r0)
orbit1.SetVelocity0(v0)
orbit1.SetIntegrationTime(t0,tf)
orbit1.SetComputeThirdBody(True)
orbit1.SetComputeSRP(False)
orbit1.SetComputeHamiltonian(False)
orbit1.SinglePropagate()
print(".")
#Get ephemeris position
ts = np.array(orbit1.getTimes())+t0
prop_pos = np.array(orbit1.getPosition()).T
ephem_pos = sp.spkpos("GRAIL-A",ts,"J2000","LT+S","MOON")[0]
err = prop_pos-ephem_pos
np.savetxt("APC prop error {:0f} to {:0f}.csv".format(t0,tf),err,delimiter=',')
# %%
