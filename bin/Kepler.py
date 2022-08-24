"""Functions for computing various properties of Keplerian Orbits and constants"""
import numpy as np

muEarth = 398600     #km^3/s^2
def rv2elm(r,v,mu,tol=1e-10):
    #  Position & Velocity Magnitudes
    R       = np.linalg.norm(r)
    V       = np.linalg.norm(v)

    #  Angular Momentum Vector
    h       = np.cross(r,v);
    H       = np.linalg.norm(h);

    #  Line of Nodes Vector
    nvec    = np.cross([0, 0, 1],h)
    n       = np.linalg.norm(nvec)

    #  Eccentricity Vector
    evec    = 1/mu*((V**2 - mu/R)*r - (np.dot(r,v))*v);
    e       = np.linalg.norm(evec);

    #  Energy
    xi      = (V**2)/2 - mu/R;

    #  Semimajor Axis (a) & Semillatus Rectum (p)
    if abs(1-e) < tol:
        a   = np.inf
        p   = (H**2)/mu
    else:
        a   = -mu/2/xi
        p   = a*(1 - e**2)

    #  Inclination
    i       = np.arccos(h[2]/H)

    # Right Ascension of Ascending Node
    Om      = np.arccos(nvec[0]/n)
    if (nvec[1] < 0):
        Om = 2*np.pi - Om

    # Argument of Perigee
    w = np.arccos((np.dot(nvec,evec))/n/e)
    if (evec[2] < 0):
        w = 2*np.pi - w

    # True Anomaly
    f = np.real(np.arccos(np.dot(evec,r)/R/e))

    if (np.dot(r,v) < 0):
        f = 2*np.pi - f

    #Mean Anomaly & Eccentric Anomaly
    E = 2*np.arctan2(np.sqrt(1-e)*np.tan(f/2),np.sqrt(1+e))
    if E < 0:
        E = 2*np.pi + E
    M = E - e*np.sin(E)
    if M < 0:
        M = 2*np.pi + M

    #  Special Cases
    #  Initialize s
    s = 0

    #  Elliptical Equatorial (ascending node undefined)
    if (i < tol) and (e >= tol):
        s = np.arccos(evec[0]/e);
        if (evec[1] < 0):
            s = 2*np.pi - s;   # Longitude of Perigee
        
    # Circular Inclined (perigee undefined)
    elif (i >= tol) and (e < tol):
        s = np.arccos(np.dot(nvec,r)/R/n);    # Argument of Latitude
        if (r[2] < 0):
            s = 2*np.pi - s;
        
    # Circular Equatorial (perigee & ascending node undefined)
    elif (i < tol) and (e < tol):
        s = np.arccos(r[0]/R);
        if (r[1] < 0):
            s = 2*np.pi - s;    # True Longitude
    return [p, a, e, i, Om, w, f, E, M, s]
def R_ECI2Orbit(Om,i,T):
    R1 = np.array([[np.cos(T),np.sin(T), 0],[-np.sin(T), np.cos(T), 0], [ 0, 0, 1]])
    R2 = np.array([[1, 0, 0],[0, np.cos(i), np.sin(i)],[0, -np.sin(i), np.cos(i)]])
    R3 = np.array([[np.cos(Om), np.sin(Om), 0],[-np.sin(Om), np.cos(Om), 0],[ 0, 0, 1]])

    R = np.matmul(R1,np.matmul(R2,R3))
    return R

def R_Orbit2ECI(Om,i,T):
    R = R_ECI2Orbit(Om,i,T).T
    return R

def elms2rv(a,e,inc,Om,w,M,mu,tol=1e-10):

    #use Netwon's method to get E
    E=M;
    g=1;
    iters=0;
    max_iter = 100000;
    while abs(g)>tol and iters<max_iter:
        g = E-e*np.sin(E)-M;
        dgde = 1-e*np.cos(E);
        E = E-g/dgde;
        iters = iters+1;
    
    #calculate true anomaly
    f = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2));

    #calculate radius and speed
    h = np.sqrt(mu*a*(1-e**2));
    #r = h**2/mu/(1+e*cos(f));
    r = a*(1-e*np.cos(E));
    T = w+f;
    #multiply r magnitude with r unit vector in ECI coordinates
    r_vec = r*np.matmul(R_Orbit2ECI(Om,inc,T),[1,0,0]);

    v_vec = mu/h*np.array([-(np.cos(Om)*(np.sin(T)+e*np.sin(w))+np.sin(Om)*(np.cos(T)+e*np.cos(w))*np.cos(inc)),
            -(np.sin(Om)*(np.sin(T)+e*np.sin(w))-np.cos(Om)*(np.cos(T)+e*np.cos(w))*np.cos(inc)),
            (np.cos(T)+e*np.cos(w))*np.sin(inc)]);
  
    return [r_vec, v_vec]

elms2rv(8000,.1,0,0,0,0,muEarth)