
from cmath import nan
import plotly.graph_objects as go
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import itertools as it

def spheres(size,pos=[0,0,0],res=100):
    #angular coordinates
    theta = np.linspace(0,2*np.pi,res)
    phi = np.linspace(0,np.pi,res)
    #cartesian coordinates
    x0 = pos[0]+size*np.outer(np.cos(theta),np.sin(phi))
    y0 = pos[1]+size*np.outer(np.sin(theta),np.sin(phi))
    z0 = pos[2]+size*np.outer(np.ones(res),np.cos(phi))

    return x0,y0,z0

def circles(size,pos=[0,0,0],res=10000):
    #angular coordinates
    theta = np.linspace(0,2*np.pi,res)
    #cartesian coordinates
    x0 = pos[0]+size*np.cos(theta)
    y0 = pos[1]+size*np.sin(theta)

    return x0,y0

def plot_sphere(size,color='green',opacity=0.3,pos=[0,0,0],res=100):
    x,y,z = spheres(size,pos,res)
    trace = go.Surface(x=x,y=y,z=z,colorscale=[[0,color],[1,color]],opacity=opacity,showscale=False)
    return trace

def plot_circle(size,color='green',opacity=0.3,pos=[0,0,0],res=100):
    x,y = circles(size,pos,res)
    cc = np.array(colors.to_rgba(color))*255
    cc[-1]=opacity
    cc=tuple(cc)
    trace = go.Scatter(x=x,y=y,line=dict(color='rgba'+str(cc)),fill='toself')
    return trace

def traceOrbit(orbit,color):
    Xs = orbit.getPosition()
    trace = go.Scatter3d(x=Xs[0][0:-2],y=Xs[1][0:-2],z=Xs[2][0:-2],mode='lines',line=dict(color=color))
    return (trace)

def traceEphem(ephem,color):
    trace = go.Scatter3d(x=ephem[:,0],y=ephem[:,1],z=ephem[:,2],mode='lines',line=dict(color=color))
    return (trace)

def plotTraces(traces,*args,**kwargs):
    if not isinstance(traces,list):
        traces = [traces]
    traceEarth = plot_sphere(1738,res=30)
    traces.append(traceEarth)
    fig = go.Figure(data=traces)
    fig.update_coloraxes(showscale=False)
    fig.update_layout(*args,**kwargs,
        scene_aspectmode='data',
        showlegend=False,
        width = 800,
        height = 800,
        margin=dict(l=100, r=100, t=100, b=100))
    return fig

def plotOrbit(orbit,color=nan,*args,**kwargs):
    if np.isnan(color):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        color=colors[0]
    trace = traceOrbit(orbit,color)
    fig = plotTraces(trace,*args,**kwargs)
    return fig

def plotEphem(ephem,color=nan,*args,**kwargs):
    if np.isnan(color):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        color=colors[0]
    trace = traceEphem(ephem,color)
    fig = plotTraces(trace,*args,**kwargs)
    return fig

def plotMultiOrbits(orbits):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = it.cycle(prop_cycle.by_key()['color'])
    for orbit in orbits:
        c = next(colors)
        traces = traceOrbit(orbit,c)
    fig = plotTraces(traces)
    return fig