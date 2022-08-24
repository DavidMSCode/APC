
import plotly.graph_objects as go
from matplotlib import colors
import numpy as np

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

def plot_sphere(size,color='blue',opacity=0.3,pos=[0,0,0],res=100):
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

def traceOrbit(orbit):
    Xs = orbit.getPosition()
    trace = go.Scatter3d(x=Xs[0],y=Xs[1],z=Xs[2],mode='lines')
    return (trace)

def plotTraces(traces):
    if not isinstance(traces,list):
        traces = [traces]
    traceEarth = plot_sphere(6378,res=30)
    traces.append(traceEarth)
    fig = go.Figure(data=traces)
    fig.update_coloraxes(showscale=False)
    fig.update_layout(
        title = "Highly eccentric orbit",
        scene_aspectmode='data',
        showlegend=False,
        width = 800,
        height = 800,
        margin=dict(l=100, r=100, t=100, b=100))
    return fig

def plotOrbit(orbit):
    trace = traceOrbit(orbit)
    fig = plotTraces(trace)
    return fig

def plotMultiOrbits(orbits):
    traces = [traceOrbit(orbit) for orbit in orbits]
    fig = plotTraces(traces)
    return fig