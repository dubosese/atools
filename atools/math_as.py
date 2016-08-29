import numpy as np
from math import cos,sin
from copy import deepcopy

def calc_com(xyz):
    xyz = np.asarray(xyz)
    return np.mean(xyz,axis=0)

def rotate(theta,dim,xyz):
    xyz = np.asarray(xyz)
    assert dim in ['x','y','z']
    com = calc_com(xyz)
    xyz -= com
    for coords in xyz:
        c = deepcopy(coords)
        if dim == 'x':
            coords[1] = c[1]*cos(theta) - c[2]*sin(theta)
            coords[2] = c[1]*sin(theta) + c[2]*cos(theta)
        elif dim == 'y':
            coords[0] = c[0]*cos(theta) + c[2]*sin(theta)
            coords[2] = -c[0]*sin(theta) + c[2]*cos(theta)
        else:
            coords[0] = c[0]*cos(theta) - c[1]*sin(theta)
            coords[1] = c[0]*sin(theta) + c[1]*cos(theta)
    xyz += com
    return xyz

def mirror(xyz):
    xyz = np.asarray(xyz)
    max = np.max(xyz[:,2])
    xyz[:,2] += (2.*(max-xyz[:,2]))
    return xyz

def func_lin(x,m,b):
    return m*x + b

def anint(x):
    if x >= 0.5:
        return 1.0
    elif x < -0.5:
        return -1.0
    else:
        return 0.0
