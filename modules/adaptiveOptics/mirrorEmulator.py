import numpy as np

class Mirror():
    def __init__(self, nPixels, nSegments, pattern=None):

        if pattern is None:
            self.pattern = np.zeros((nPixels,nPixels))
        else:
            self.pattern = pattern

        self.nPixels = nPixels
        self.nSegments = nSegments

        self.segOffsets = np.zeros((nSegments, nSegments))
        self.segTilts = np.zeros((nSegments, nSegments))

        self.borders = np.linspace(0,self.nPixels,num=nSegments+1)

        self.pixInSeg = self.borders[1]-self.borders[0]

        self.setGeometry(nPixels/2,nPixels/2,nPixels)

    def setGeometry(self, cx, cy, d):
        '''
        Sets geometry on pixellated grid
        '''
        self._geometry = pupil.Geometry((self.nPixels,
                                         self.nPixels),
                                        cx,cy,d)
        


    def whereSegment(self, seg):
        temp = np.zeros_like(self.pattern)
        unraveled = np.unravel_index(seg, [self.nSegments, self.nSegments])
        xStart = self.borders[unraveled[0]]
        xStop = self.borders[unraveled[0]+1]
        yStart = self.borders[unraveled[1]]
        yStop = self.borders[unraveled[1]+1]
        temp[xStart:xStop,yStart:yStop] = -1
        return np.where(temp==-1)

    def addOffset(self, seg, value):
        unraveled = np.unravel_index(seg, [self.nSegments, self.nSegments])
        self.segOffsets[unraveled[0],unraveled[1]] += value
        w = self.whereSegment(seg)
        self.pattern[w] += seg

    def findSegOffsets(self):
        for i in range(self.nSegments*self.nSegments):
            w = self.whereSegment(i)
            av = np.mean(self.pattern[w])
            unraveled = np.unravel_index(i, [self.nSegments, self.nSegments])
            self.segOffsets[unraveled[0],unraveled[1]]=av

    def findTilt(self, tilt):
        temp = tilt * np.arange(self.pixInSeg)
        toAdd = temp.repeat(self.pixInSeg).reshape(self.pixInSeg,self.pixInSeg)
        return toAdd

    def addXTilt(self, seg, tilt):
        w = self.whereSegment(seg)
        toAdd = self.findTilt(tilt)
        print "toAdd shape:", toAdd.shape
        print "pattern segment shape:", self.pattern[w].shape
        self.pattern[w] += toAdd
        
    def addYTilt(self, seg, tilt):
        w = self.whereSegment(seg)
        toAdd = self.findTilt(tilt)
        self.pattern[w] += toAdd.swapaxes(0,1)

    def calcZernikes(self, zernModes):
        radius = self.nPixels/2
        zernPattern = libtim.zern.calc_zernike(zernModes, radius, mask=False)
        return zernPattern

    def clear(self):
        self.pattern = np.zeros((self.nPixels,self.nPixels))
        self.segOffsets = np.zeros((self.nSegments, self.nSegments))

    def addZernikes(self, zernModes, clear_first = False):
        if clear_first:
            self.clear()
        z = self.calcZernikes(self, zernModes)
        self.pattern += z
        self.findSegOffsets()

    def pokeSegment(self, seg, amount):
        pass

    def outputSegs(self, filename):
        segs = np.append(self.segOffsets.astype(np.int16).flatten(),np.zeros((16),dtype=np.int16))
        np.savetxt(filename, segs, fmt='%i', delimiter='\n')
        return segs
        
        
if __name__ == "__main__":

    pixels = 512
    segments = 12
    toSaveDir = "Z:\\Ryan\\BMC_Mirror\\seg\\"
    m = Mirror(pixels, segments)
    allFiles = []
    for i in range(0,segments*segments):
        m.addOffset(i,1000)
        fname = toSaveDir+"poked"+str(i).zfill(3)+".txt"
        m.outputSegs(fname)
        allFiles.append(fname)
        m.clear()
    np.savetxt(toSaveDir+"allSegFiles.txt", np.array(allFiles), fmt='%s', delimiter='\n')
