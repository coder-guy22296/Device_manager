import numpy as np
import scipy
from scipy import optimize
from math import hypot, atan2, pi
import os

def polar(x,y,deg=0):
    if deg:
        return hypot(x,y), 180.0*atan2(y,x) / pi
    else:
        return hypot(x,y), atan2(y,x)


def polar_numpy(x,y,deg=0):
    if deg:
        return np.sqrt(x**2 + y**2), 180.0*np.arctan2(y,x) / pi
    else:
        return np.sqrt(x**2 + y**2), 180.0*np.arctan2(y,x)

    
def fitPlane(pf, NAmask):
    w = np.where(NAmask>0)
    nx,ny = pf.shape
    xs,ys = np.meshgrid(np.arange(-1*nx/2,nx/2),np.arange(-1*ny/2,ny/2))
    points = np.vstack((xs[w[0],w[1]].flatten(), ys[w[0],w[1]].flatten(), pf[w[0],w[1]].flatten()))
    ctr = points.mean(axis=1)
    x = points - ctr[:,None]
    M = np.dot(x, x.T)
    solution = np.linalg.svd(M)[0][:,-1]
    fittedPlane = (-1*solution[0]*(xs-ctr[0]) - solution[1]*(ys-ctr[1])) / solution[2]
    return (pf-fittedPlane)*NAmask

def fitDefocus(pf, NAmask):
    nx,ny = pf.shape
    xs,ys = np.meshgrid(np.arange(-1*nx/2,nx/2),np.arange(-1*ny/2,ny/2))
    r,t = polar_numpy(xs,ys)
    p0 = [-0.001, (pf*NAmask).mean()]
    p1, success = optimize.leastsq(errorInDefocus, p0, args=(pf*NAmask, r, NAmask))
    #print "Fitting defocus result: ", (p1, success)
    if success != 1:
        print "Something not right with defocus fit..."
    fitResult = (r**2*p1[0] + p1[1])*NAmask
    return pf - fitResult

def errorInDefocus(params, pf, r, NAmask):
    guess = (r**2*params[0] + params[1])*NAmask
    error = guess - pf
    return error.flatten()

