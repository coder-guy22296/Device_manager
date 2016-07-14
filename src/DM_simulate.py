# simulate the deformable mirror 
import numpy as np 
import libtim.zern 
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from scipy import optimize, integrate
import pupil2device as pupil
from numpy.lib.scimath import sqrt as _msqrt

import time


class DM_simulate(object):
    def __init__(self, nseg = 12, nPixels = 256, pattern=None):
        self.nSegments = nseg
        self.nPixels = nPixels
        self.segOffsets = np.zeros((self.nSegments, self.nSegments))
        
        if pattern is None:
            self.pattern = np.zeros((nPixels,nPixels))
        else: 
            self.pattern = pattern
            
        self.initGeo(nPixels)    
    
            
    def discretize(self,weighted = None):
        pass 
    
    
    def zernike_single(self, mode, amp, rad = None):
        modes = np.zeros((mode))
        modes[mode-1]=amp
        if rad is None:
            radius = self.nPixels/2
        else:
            radius = rad
            
            
        zernike = libtim.zern.calc_zernike(modes, radius, mask=True,zern_data = {})
        self.pattern = zernike
        return zernike 

    
    
    def initGeo(self, npixels = 256):
        self.nPixels = npixels
        self.borders = np.linspace(0,self.nPixels,num=self.nSegments+1).astype(int)
        self.pixInSeg = self.borders[1]-self.borders[0]
        return self.borders


    def findSegOffsets(self, pat):
        # This is where segment values are updated
             
        for i in range(self.nSegments*self.nSegments):
            w = self.whereSegment(i)
            av = np.mean(pat[w])
#             print(av)
            unraveled = np.unravel_index(i, [self.nSegments, self.nSegments])
            self.segOffsets[unraveled[0],unraveled[1]]=av
            
        segOffsets = np.copy(self.segOffsets)    
        return segOffsets
        
            
    def whereSegment(self, seg):
        temp = np.zeros_like(self.pattern)
        unraveled = np.unravel_index(seg, [self.nSegments, self.nSegments])
        xStart = self.borders[unraveled[0]]
        xStop = self.borders[unraveled[0]+1]
        yStart = self.borders[unraveled[1]]
        yStop = self.borders[unraveled[1]+1]
#         print(xStart, xStop, yStart, yStop)
        temp[xStart:xStop,yStart:yStop] = -1
        return np.where(temp==-1)

    def findSeg(self):
        for ii in np.arange(self.nSegments):
            for jj in np.arange(self.nSegments):
                xStart = self.borders[ii]
                xEnd = self.borders[ii+1]
                yStart = self.borders[jj]
                yEnd = self.borders[jj+1]
                
                av = np.mean(self.pattern[xStart:xEnd, yStart:yEnd])
                
#                 w = self.whereSeg(ii, jj)
#                 av = np.mean(self.pattern[w])
                self.segOffsets[ii,jj] = av
                
        segOffsets = np.copy(self.segOffsets)
        return segOffsets
        
        
        
    def whereSeg(self, ii, jj):
        xStart = self.borders[ii]
        xEnd = self.borders[ii+1]
        yStart = self.borders[jj]
        yEnd = self.borders[jj+1]
        
        xRange = np.arange(xStart,xEnd)
        yRange = np.arange(yStart,yEnd)
        
        MY, MX = np.meshgrid(yRange, xRange) 
        block_ind = (MX.ravel(), MY.ravel())
        return block_ind
        
        
        
        

def main():
    ndeg = 22
    inner_prod = np.zeros([ndeg,ndeg])
    DM = DM_simulate(nseg = 256, nPixels = 256)
    for N1 in np.arange(ndeg):
        zen_1 = DM.zernike_single(N1+1, 1.0)
#         seg_1 = DM.findSegOffsets(zen_1)
        seg_1 = DM.findSeg()
        seg_1 = seg_1.ravel()
        for N2 in np.arange(N1+1):
            zen_2 = DM.zernike_single(N2+1, 1.0)
            seg_2=DM.findSeg()
#             seg_2 = DM.findSegOffsets(zen_2)
            seg_2 = seg_2.ravel()
            inner_prod[N1, N2] = np.inner(seg_1, seg_2) / len(seg_1)
            inner_prod[N2, N1] = np.inner(seg_1, seg_2) / len(seg_1)#     pylab.imshow(seg_1)
            print(N1,N2)
            
    im = plt.imshow(inner_prod,interpolation='none',cmap = 'RdBu')
    plt.colorbar(im)
    plt.savefig('seg_256_new')
#     plt.imshow()
#     print(inner_prod)


if __name__ == "__main__":
    stime = time.time()
    main()
    etime = time.time()-stime
    print(etime)
# print(zen_1.shape)
# pylab.imshow(zen_1)
# pylab.show()
