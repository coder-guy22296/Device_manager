import inLib
import time
import numpy as np
import libtim
import libtim.zern
from modules.adaptiveOptics import pupil_forInControl as pupil
import subprocess
import scipy
from scipy.ndimage import rotate
from Cython.Compiler.ExprNodes import constant_value_not_set

class Control(inLib.Device):
    def __init__(self, settings):

        self.pixels = settings['pixels']
        self.segments = settings['segments']
        self.executable = settings['executable']

        self.mirror = Mirror(self.pixels, self.segments)
        self.pad_mirror = Mirror(self.pixels, self.segments)
        self.dummy_mirror = Mirror(256,12)

        self._geometry = self.mirror.geometry
        self.tempfilename = 'mirrorSegs.txt'
        self.multiplier = 1.0
        self.preMultiplier = 80
        
        self.group = []
        self.padding = False

        self.proc = None

        #below is added by Dan 
        # have a stack of zernike modes
        self.z_max = 25
        self.zernike = np.zeros(self.z_max) # Initialize self.zernike 
        self.pool = Modulation_pool(self.z_max)
        

    def getGeometry(self):
        return self._geometry

    def getGeoParams(self):
        '''Returns npixels, cx, cy'''
        return self.mirror.getGeo()
    
    def highlight_dummy_mirror_segs(self, segs):
        self.group = segs
        for i in range(len(segs)):
            self.dummy_mirror.addOffset(segs[i], 1)
        return self.dummy_mirror.pattern


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


    def setMultiplier(self,mult, pt = None):
        if pt is None:
            self.multiplier = mult
        else:
            self.multi_list[pt] = mult
            

    def setPreMultiplier(self,mult):
        self.preMultiplier = mult    


    def returnSegments(self):
        return self.mirror.returnSegs()
        #return self.mirror.segOffsets

    def returnPattern(self):
        return self.mirror.pattern
        #But how is this pattern update from addOther?
        

    def clear(self):
        self.mirror.clear()

    
    def calcZernike_single(self, mode, amp, radius=None, useMask=True):
        # this is for single Zernike
        if radius is None:
            radius = self.mirror.nPixels/2
        modes = np.zeros((mode))
        modes[mode-1]=amp
        self.mirror.pattern = libtim.zern.calc_zernike(modes, radius, mask=useMask,
                                                zern_data = {})
        return self.mirror.pattern
    
    def calcZernike_multi(self, amps, radius = None, useMask = True):
        if radius is None:
            radius = self.mirror.nPixels/2
        self.mirror.pattern = libtim.zern.calc_zernike(amps,radius, mask = useMask, zern_data = {})
        
        return self.mirror.pattern
            
    # disabled on 07/20
#     def addZernike(self, zernike_pattern=None):
#         if zernike_pattern is not None:
#             zern = zernike_pattern
#         else:
#             if self.zernike is None:
#                 return 0
#             else:
#                 zern = self.zernike
#         if self.zernike is not None:
#             self.mirror.addToPattern(zern * self.preMultiplier) # Here premultiplyer already counts for part of the scaling factor.
# #             zernike_pattern            
#             
#         return self.returnPattern()

  
    def advancePatternWithPipe(self):
        if self.proc is not None:
            self.proc.stdin.write("\n")

    def resetPipe(self):
        self.proc = None


    def applyToMirror(self):
        #First save mirror
        self.mirror.outputSegs(self.tempfilename)
        #Wait to make sure file exists
        time.sleep(0.5)
        #Run executable
        subprocess.call([self.executable, self.tempfilename, str(self.multiplier),
                         "1", "-1"], shell=True)
                         
                         
    def addOther(self, MOD):
        # Update by Dan on 07/14. 
        pass
    
    # ------------------- everything of zernike mode modulation and pool operation---------------------
    
    
    def push_to_zernike(self, nmode, amp):
        # push a single mode into self.zernike
        self.zernike[nmode-1] = amp 
        print("added mode:", nmode)
        
    
    def push_to_pool(self, zm, multi):
        # added on 07/20: nmode is the zm, amp plays the multiplier's role
        # then clear self.zernike for the next group of modulation
        self.pool.append_mod(zm, multi)
        self.zernike = np.zeros(self.z_max)
    
    
    def setMod_status(self,index,state):
        # this feels so awkward to handle. But let it be.
        self.pool.z_active[index] = state
        
        
    def mod_from_pool(self):
        amps = self.pool.synth_mod()
        self.mirror.pattern = self.calcZernike_multi(amps)
        self.mirror.findSeg()
        # Up to here, the pattern is not added to the mirror yet.
#         self.applyToMirror()
        
        
    
    
class Modulation_pool(object):
    def __init__(self,z_max):
        self.multi_list = []
        self.z_max = z_max
        self.zen_list = []
        self.z_store = 0
        self.z_active = []
    
    def append_mod(self, zm, multi=1):
        npad = self.z_max - len(zm)
        if npad>0:
            zmodes = np.lib.pad(zm,(0, npad), 'constant', constant_value = 0.)
        else:
            zmodes = zm[:self.z_max] 
        self.zen_list.append(zmodes)
        self.multi_list.append(multi) 
        self.z_store +=1
        self.z_active.append(True)
        
    def delete_mod(self, ndel = -1):
        # by default: delete the last element 
        del(self.zen_list[ndel])
        del(self.multi_list[ndel])
        del(self.z_active[ndel])
        self.z_store -=1
        
        
    def synth_mod(self):
        # find the Active modulations
        n_active = [i for i, s in enumerate(self.z_active) if s] # select all the True modes
        z_coeff = np.zeros_like(self.zen_list[0], dtype='float')
        for nmod in n_active:
            z_coeff += self.multi_list[nm]*self.zen_list[nm]
            
        return z_coeff 
        """ 
        This is the final zernike modes scaled by all the multipliers. When applying to the mirror, there's no
        need to rescale it with
        """    
            
        

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

            
    def findSeg(self):
        # This is a much simpler version than findSegOffsets
        for ii in np.arange(self.nSegments):
            for jj in np.arange(self.nSegments):
                xStart = self.borders[ii]
                xEnd = self.borders[ii+1]
                yStart = self.borders[jj]
                yEnd = self.borders[jj+1]
                
                av = np.mean(self.pattern[xStart:xEnd, yStart:yEnd])
                self.segOffsets[ii,jj] = av
                


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


    def clear(self):
        self.pattern = np.zeros((self.nPixels,self.nPixels))
        self.segOffsets = np.zeros((self.nSegments, self.nSegments))


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
