from numpy import *
from scipy import optimize

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def gaussian1d(height, center_x, width_x):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    return lambda x: height*exp(-1*(((center_x-x)/width_x)**2))

def moments1d(data):
    """Returns (height, x, width_x)
    the gaussian parameters of a 1D distribution by calculating its
    moments """
    total = data.sum()
    X = arange(len(data))
    x = (X*data).sum()/total
    width_x = sqrt(abs((X-x)**2*data).sum()/data.sum())
    height = data.max()
    return height, x,width_x

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fitgaussian1d(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments1d(data)
    errorfunction = lambda p: ravel(gaussian1d(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p
