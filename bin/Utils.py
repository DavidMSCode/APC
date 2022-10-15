import numpy as np
from Kepler import *

def out2breaks(output):
    Ts = output.getTimes()
    Xs = output.getPosition()
    Vs = output.getVelocity()
    #Segment orbits into individual revolutions of the Earth
    prev_y=0.0
    orbit_breaks =[0]
    for i,y in enumerate(Xs[1]):
        if prev_y<0 and y>0:
            orbit_breaks.append(i)
        prev_y = y
    orbits_x = []
    orbits_y = []
    orbits_t = []
    for i,start in enumerate(orbit_breaks[0:-1]):
        if not start==0:
            start=start-1
        end = orbit_breaks[i+1]+1
        orbits_x.append(Xs[0][start:end])
        orbits_y.append(Xs[1][start:end])
        orbits_t.append(Ts[start:end])

    #find orbital elements at each time step
    elms = np.zeros((len(orbit_breaks),10))
    ts = np.zeros(len(orbit_breaks))
    for i,ind in enumerate(orbit_breaks):
        x = np.array([Xs[0][ind],Xs[1][ind],Xs[2][ind]])
        v = np.array([Vs[0][ind],Vs[1][ind],Vs[2][ind]])
        elms[i,:] = rv2elm(x,v,muEarth)
        #[p, a, e, i, Om, w, f, E, M, s]
    
    return [orbits_x, orbits_y, orbits_t, elms] 