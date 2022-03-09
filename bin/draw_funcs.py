
import plotly.graph_objects as go
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
    trace = go.Scatter(x=x,y=y,line=dict(color=color),fill='toself')
    return trace