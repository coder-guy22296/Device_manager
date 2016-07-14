import inLib
import time
import numpy as np
import libtim
import libtim.zern
from modules.adaptiveOptics import pupil_forInControl as pupil
import subprocess
import scipy
from scipy.ndimage import rotate

class Control(inLib.Device):
    def __init__(self, settings):

        self.pixels = settings['pixels']
        self.segments = settings['segments']
        self.executable = settings['executable']

        self.mirror = Mirror(self.pixels, self.segments)

        self.pad_mirror = Mirror(self.pixels, self.segments)

        self.dummy_mirror = Mirror(256,12)
        self.setup_dummy_mirror()

        self._geometry = self.mirror.geometry

        self.tempfilename = 'mirrorSegs.txt'
        self.multiplier = 1.0

        self.preMultiplier = 80

        self.zernike = None

        self._other = _OtherModulations()

        self.group = []

        self.padding = False

        self.proc = None
        self.numZernsToVary = 100

    def getGeometry(self):
        return self._geometry

    def getGeoParams(self):
        '''Returns npixels, cx, cy'''
        return self.mirror.getGeo()

    def setup_dummy_mirror(self):
        pass

    def highlight_dummy_mirror_segs(self, segs):
        self.group = segs
        for i in range(len(segs)):
            self.dummy_mirror.addOffset(segs[i], 1)
        return self.dummy_mirror.pattern

    def clearDummyMirror(self):
        self.dummy_mirror.clear()

    def returnDummyMirror(self):
        return self.dummy_mirror.pattern, self.dummy_mirror.segOffsets
        

    def loadPattern(self, pattern_filename, mult=1.0):
        self.mirror.pattern = np.load(str(pattern_filename)) * mult
        return self.mirror.pattern

    def loadSegments(self, filename):
        msegs = np.loadtxt(filename)
        #print "msegs: ", msegs
        self.mirror.inputMirrorSegs(msegs)
        
    def patternRot90(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.rot90(self.mirror.pattern)
        return self.mirror.pattern

    def patternFlipLR(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.fliplr(self.mirror.pattern)
        return self.mirror.pattern

    def patternFlipUD(self):
        if self.mirror.pattern is not None:
            self.mirror.pattern = np.flipud(self.mirror.pattern)
        return self.mirror.pattern

    def patternRotate(self, deg):
        if self.mirror.pattern is not None:
            self.mirror.pattern = rotate(self.mirror.pattern, deg)
        return self.mirror.pattern

    def pokeGroup(self, group, offset, quiet=False):
        for s in group:
            self.pokeSegment(s, offset)
            
    def pokeSegment(self, seg, value, pokeAll = False):
        '''
        Add value to a given segment or to all segments
        '''
        if pokeAll:
            self.mirror.addOffset(-1,value)
        else:
            self.mirror.addOffset(seg,value)

    def reconfigGeo(self, cx, cy, npixels):
        self.mirror.initGeo(npixels)
        half_px = npixels/2
        if self.mirror.pattern is not None:
            if self.mirror.pattern.shape[0] > npixels:
                pattern = self.mirror.pattern[cx-half_px:cx+half_px,
                                              cy-half_px:cy+half_px]
                if (pattern.shape[0] == npixels) and (pattern.shape[1]==npixels):
                    self.mirror.pattern = pattern.copy()
                else:
                    print "Problem cropping pattern for DM."
            if self.mirror.pattern.shape[0] < npixels:
                pattern = np.zeros((npixels,npixels))
                origShape = self.mirror.pattern.shape[0]
                pattern[cx-origShape/2:cx+origShape/2, cy-origShape/2:cy+origShape/2] = self.mirror.pattern
                self.mirror.pattern = pattern.copy()
        self._geometry = self.mirror.geometry
        return self.returnPattern()

    def padZeros(self, border, always=False):
        '''
        Pad pattern with zeros
        '''
        currentPixels = self.mirror.pattern.shape[0]
        newPattern = np.zeros((currentPixels+2*border, currentPixels+2*border))
        newPattern[border:-1*border,border:-1*border] = self.mirror.pattern
        self.mirror.initGeo(newPattern.shape[0])
        self.mirror.setPattern(newPattern)
        return self.returnPattern()

    def findSegments(self):
        self.mirror.findSegOffsets()

    def setMultiplier(self,mult):
        self.multiplier = mult

    def setPreMultiplier(self,mult):
        self.preMultiplier = mult

    def getSegments(self):
        return self.mirror.returnSegs()
        #return self.mirror.segOffsets

    def returnSegments(self):
        return self.mirror.returnSegs()
        #return self.mirror.segOffsets

    def returnPattern(self):
        return self.mirror.pattern

    def clear(self):
        self.mirror.clear()

    def setZernMode(self, mode):
        self.zernMode = mode
    
    def calcZernike(self, mode, amp, radius=None, useMask=True):
        if radius is None:
            radius = self.mirror.nPixels/2
        modes = np.zeros((mode))
        modes[mode-1]=amp
        self.zernike = libtim.zern.calc_zernike(modes, radius, mask=useMask,
                                                zern_data = {})
        return self.zernike

    def addZernike(self, zernike_pattern=None):
        if zernike_pattern is not None:
            zern = zernike_pattern
        else:
            if self.zernike is None:
                return 0
            else:
                zern = self.zernike
        if self.zernike is not None:
            p=self.mirror.addToPattern(zern * self.preMultiplier)
        return self.returnPattern()

    def getNumberOfZernToVary(self):
        return self.numZernsToVary

    def varyMultiplierCurrent(self, minMult, maxMult, num, wTime, externallyCalled = False):
        print("The number of multipliers:", num)
        print("mininum:", minMult)
        print("maximum:",maxMult)
        self.numZernsToVary = num
        mults = np.linspace(minMult,maxMult,num)
        filenms = []
        baseline = self.returnPattern()
        for i in xrange(num):
            self.clear()
            newPattern = baseline * mults[i]
            self.mirror.setPattern(newPattern)
            self.findSegments()
            filenms.append("segFile_mult_%.2i.txt" % i)
            self.mirror.outputSegs(filenms[i])
        files_file = "allMultFiles_Max%.2i.txt" % maxMult
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print "Finished creating files for varying multiplier for current pattern..."

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr]

        if self.proc is not None:
            print "Polling proc: ", self.proc.poll()
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print "Going to next..."
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print "proc stdout: ", output
            return 1
        

    def varyZernAmp(self, mode, maxAmp, minAmp, num, wTime, radius=None, useMask=True, clearfirst=True,
                    externallyCalled = False):
        '''
        Calls external C++ program that applies *num* different patterns of
        Zernike mode *mode* to the mirror.

        :Parameters:
            *mode*: int
                Zernike mode
            *maxAmp*: float
            *minAmp*: float
            *num*: int
            *wTime*: float
                Time to wait in milliseconds before new pattern applied to mirror
            *useMask*: boolean
                Optional
            *clearFirst*: boolean
                Optional
        '''
        #mode = self.zernMode
        self.numZernsToVary = num
        amps = np.linspace(minAmp,maxAmp,num)
        filenms = []
        baseline = self.returnPattern()
        for i in xrange(num):
            if clearfirst:
                self.clear() #clears mirror pattern and segments
            else:
                self.clear()
                self.mirror.setPattern(baseline)
            zern = self.calcZernike(mode, amps[i], radius=radius, useMask=useMask)
            self.addZernike(zernike_pattern=zern)
            self.findSegments()
            filenms.append("segFile_mode%.3i_amp%.2i.txt" % (mode,i))
            self.mirror.outputSegs(filenms[i])
        files_file = "allSegFiles_%.3i.txt" % mode
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print "Finished creating files for varying Zernike..."

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr] 

        #subprocess.call([self.executable, files_file, str(self.multiplier),
        #                 numStr, wTimeStr], shell=True)

        if self.proc is not None:
            print "Polling proc: ", self.proc.poll()
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print "Going to next..."
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print "proc stdout: ", output
            return 1

    def varyZernRadii(self, mode, amp, maxR, minR, num, wTime, radius=None, useMask=True, clearfirst=True,
                      externallyCalled = False):
        '''
        Calls external C++ program that applies *num* different patterns of
        Zernike mode *mode* to the mirror.

        :Parameters:
            *mode*: int
                Zernike mode
            *amp*: float
                Zernike amplitude
            *maxR*: float
            *minR*: float
            *num*: int
            *wTime*: float
                Time to wait in milliseconds before new pattern applied to mirror
            *useMask*: boolean
                Optional
            *clearFirst*: boolean
                Optional
        '''
        #mode = self.zernMode
        self.numZernsToVary = num
        rads = np.linspace(minR,maxR,num,dtype=np.uint16)
        filenms = []
        baseline = self.returnPattern()
        for i in xrange(num):
            if clearfirst:
                self.clear() #clears mirror pattern and segments
            else:
                self.clear()
                self.mirror.setPattern(baseline)
            zern = self.calcZernike(mode, amp, radius=rads[i], useMask=useMask)
            self.addZernike(zernike_pattern=zern)
            self.findSegments()
            filenms.append("segFile_mode%.3i_rad%.2i.txt" % (mode,i))
            self.mirror.outputSegs(filenms[i])
        files_file = "allSegFiles_%.3i.txt" % mode
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        time.sleep(0.1)

        print "Finished creating files for varying Zernike radii..."

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr] 

        #subprocess.call([self.executable, files_file, str(self.multiplier),
        #                 numStr, wTimeStr], shell=True)

        if self.proc is not None:
            print "Polling proc: ", self.proc.poll()
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print "Going to next..."
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print "proc stdout: ", output
            return 1
        

    def advancePatternWithPipe(self):
        if self.proc is not None:
            self.proc.stdin.write("\n")

    def resetPipe(self):
        self.proc = None

    def varyArbitrary(self, filename, wTime, externallyCalled=False):
        wTimeStr = "%i" % wTime
        files_file = str(filename)
        num = len(np.loadtxt(files_file, dtype='S'))
        print "varyArbitrary: Number of files to load is ", num
        numStr = "%i" % num
        self.numZernsToVary = num
        args = [self.executable, files_file, str(self.multiplier), numStr, wTimeStr]

        if self.proc is not None:
            print "Polling proc: ", self.proc.poll()
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None

        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print "Going to next..."
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print "proc stdout: ", output
            return 1
        
            
    def varyGroupOffset(self, maxAmp, minAmp, num, wTime, group=None, useMask=True,
                        clearfirst=True, externallyCalled = False):
        amps = np.linspace(minAmp,maxAmp,num)
        filenms = []
        self.numZernsToVary = num
        baseline = self.returnSegments()
        if group is None:
            group = self.group
        for i in range(num):
            #Need to change clearfirst
            if clearfirst:
                self.clear() #clears mirror pattern and segments
            else:
                self.clear()
                self.mirror.useSegements(baseline)
            self.pokeGroup(group, amps[i])
            #self.findSegments() #don't need. we're directly adding segments
            filenms.append("segFile_GroupVary_amp%.2i.txt" % i)
            self.mirror.outputSegs(filenms[i])
        files_file = "allSegFiles_group.txt"
        np.savetxt(files_file, np.array(filenms), fmt='%s', delimiter='\n')

        print "Finished creating files for varying group offsets..."
        time.sleep(0.5)

        wTimeStr = "%i" % wTime
        numStr = "%i" % num

        args = [self.executable, files_file, str(self.multiplier),numStr, wTimeStr]
        
        if self.proc is not None:
            print "Polling proc: ", self.proc.poll()
            if self.proc.poll() is None:
                self.proc.terminate()
                self.proc.communicate()
                self.proc = None
            
        self.proc = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        if externallyCalled:
            return 0
        else:
            if wTime<0:
                for i in range(num):
                    time.sleep(1)
                    print "Going to next..."
                    self.advancePatternWithPipe()
            output = self.proc.stdout.read()
            print "proc stdout: ", output
            return 1


    def addOther(self, data, on_geometry_changed=None, args=()):
        data_new = np.rot90(data-data.min())
        print "Data_New max: ", data_new.max()
        data_new = self.preMultiplier * data_new #/ data_new.max()
        index = self._other.counter
        self._other.others[index] = _OtherModulation(data_new-data_new.mean(),on_geometry_changed, args)
        self._other.counter += 1
        #if self._active and self._other.active:
        #    self._sum_send()
        return index

    def setOtherActive(self, index, state):
        self._other.others[index].active = state
        if state:
            p=self.mirror.addToPattern(self._other.others[index].data)
        return self.returnPattern()

    def applyToMirror(self, wtime=-1):
        #First save mirror
        self.mirror.outputSegs(self.tempfilename)

        #Wait to make sure file exists
        time.sleep(0.5)

        wTimeStr = str(wtime)

        #Run executable
        subprocess.call([self.executable, self.tempfilename, str(self.multiplier),
                         "1", wTimeStr], shell=True)
        #subprocess.call([self.executable, self.tempfilename, str(self.multiplier),
         #                "1", "30000"], shell=True)
        

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

        #self.setPupilGeometry(nPixels/2,nPixels/2,nPixels)

        self.initGeo(self.nPixels)

    def initGeo(self, npixels):
        self.nPixels = npixels
        self.borders = np.linspace(0,self.nPixels,num=self.nSegments+1)
        self.pixInSeg = self.borders[1]-self.borders[0]
        self.setPupilGeometry(npixels/2, npixels/2, npixels)
        

    def setPupilGeometry(self, cx, cy, d):
        '''
        Sets geometry on pixellated grid
        '''
        self.geometry = pupil.Geometry((self.nPixels,
                                         self.nPixels),
                                        cx,cy,d)

    def getGeo(self):
        return [self.geometry.nx, self.geometry.cx, self.geometry.cy]


    def setPattern(self, data):
        self.pattern = data.copy()
        return self.pattern

    def addToPattern(self, data):
        if self.pattern.shape == data.shape:
            self.pattern += data
        elif data.shape < self.pattern.shape:
            diffx = self.pattern.shape[0] - data.shape[0]
            diffy = self.pattern.shape[1] - data.shape[1]
            if diffx != diffy:
                print "Mirror.addToPattern: Something's not square..."
            else:
                border = diffx/2
                self.pattern[border:-1*border,border:-1*border] += data
        else:
            print "Mirror.addToPattern: Shape mismatch..."
            print "Mirror.addToPattern: Data to add has shape: ", data.shape
            print "Mirror.addToPattern: Pattern  has shape: ", self.pattern.shape
        return self.pattern
        

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
        #print "adding %i to segment %i." % (value,seg)
        if seg<0:
            self.segOffsets += value
            self.pattern += value
        else:
            unraveled = np.unravel_index(seg, [self.nSegments, self.nSegments])
            self.segOffsets[unraveled[0],unraveled[1]] += value
            w = self.whereSegment(seg)
            self.pattern[w] += value

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

    '''
    def calcZernikes(self, zernModes):
        radius = self.nPixels/2
        zernPattern = libtim.zern.calc_zernike(zernModes, radius, mask=False)
        return zernPattern
    '''  

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
        allSegments = self.segOffsets.astype(np.int16).flatten()
        forMirror = np.zeros((160),dtype=np.int16)
        forMirror[0:10] = allSegments[1:11]
        forMirror[10:130] = allSegments[12:132]
        forMirror[130:140] = allSegments[133:143]
        segs = np.append(forMirror, np.zeros((16),dtype=np.int16))
        #segs = np.append(self.segOffsets.astype(np.int16).flatten(),np.zeros((16),dtype=np.int16))
        np.savetxt(filename, segs, fmt='%i', delimiter='\n')
        return segs

    def inputMirrorSegs(self, msegs):
        newSegs = np.zeros((144))
        newSegs[1:11] = msegs[0:10]
        newSegs[12:132] = msegs[10:130]
        newSegs[133:143] = msegs[130:140]
        self.segOffsets = newSegs.reshape(12,12)
        #print "segOffsets: ", self.segOffsets

    def returnSegs(self):
        return self.segOffsets

    def useSegements(self,segs):
        self.segOffsets = segs.copy()

class _OtherModulations:
    def __init__(self):
        self.active = True
        self.others = {}
        self.counter = 0

class _OtherModulation:

    def __init__(self, data, on_geometry_changed, args):

        self.active = True
        self.data = data
        self.on_geometry_changed = on_geometry_changed
        self.args = args
        
        
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
