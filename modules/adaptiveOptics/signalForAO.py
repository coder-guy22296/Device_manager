#!/usr/bin/python
#
## @file
#
# 
# RJM 2/15


import numpy as np
import os.path
import time
import pyfftw
import sys

def fftImage(image, use_pyfftw=True):
    if use_pyfftw:
        pyfftw.interfaces.cache.enable()
        im = pyfftw.n_byte_align(image, 16)
        return pyfftw.interfaces.numpy_fft.fftshift(pyfftw.interfaces.numpy_fft.fft2(image))
    return np.fft.fftshift(np.fft.fft2(image))

def window(ndim, pixelSize, diffLimit):
    x,y = np.meshgrid(np.arange(-1*ndim/2,ndim/2),
                      np.arange(-1*ndim/2,ndim/2))
    r = np.sqrt(x**2 + y**2)
    freqs = r/(ndim*pixelSize)
    res = 1./diffLimit
    return freqs,freqs<res

def secondMoment(image, pixelSize, diffLimit):
    ndim = image.shape[0]
    freqs, win = window(ndim,pixelSize,diffLimit)
    fI = np.abs(fftImage(image))
    numerator = np.sum(win*fI*(freqs**2))
    return numerator/np.sum(fI)

def secondMomentOnStack(images, pixelSize, diffLimit):
    if images is not None:
        sm = np.zeros((images.shape[0]))
        for i in range(0,len(sm)):
            sm[i] = secondMoment(images[i],pixelSize,diffLimit)
        return sm
    else:
        print "got None for images"
