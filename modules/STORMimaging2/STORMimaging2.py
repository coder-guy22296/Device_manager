import inLib
import time
import numpy as np
from scipy.fftpack import fft2,ifft2,fftshift
import libs.imagewriters as writer
import os, glob


class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing STORM imaging..."
        inLib.Module.__init__(self, control, settings)
        print "STORM imaging initialized."

        self.initSettings = settings
        self.initSettings['camera'] = {}
        self.initSettings['camera']['settings'] = self._control.camera._initSettings
        self.initSettings['settings_filename'] = "Initial settings"
        

        #self._control.piezo ==> controls piezo
        #self._control.shutters ==> controls NIDAQ
        #self._control.camera ==> Hamamatsu
        #self._control.stage ==> Marhauser

        self.savedLocations = []
        self.savedSettings = [self.initSettings]

        self.daxfile = 0
        self.filming = 0
        self.zerothFrame = 0
        self.fixedLength = False

        self.openShutterState = False


    def reInit(self):
        self.all_data = []
        self.stage_zs = []
        self.times = []
        self.means = []

    def getXYZPosition(self):
        x,y,z = self._control.stage.position()
        z = self._control.piezo.getPosition(3)
        return x,y,z

    def getXYFocus(self):
        x,y,z = self.getXYZPosition()

    def moveAndSnap(self, dx, dy):
        imnum = self._control.camera.takeSnap()
        self._control.stage.goRelative(dx,dy)
        return imnum

    def rememberLocation(self):
        x,y,z = self.getXYZPosition()
        if len(self.savedLocations)<4:
            self.savedLocations.append([x,y])
        return len(self.savedLocations)

    def clearLocations(self):
        self.savedLocations = []

    def goToLocation(self, num):
        if num<=len(self.savedLocations):
            x,y = self.savedLocations[num-1]
            print "Moving to %.1f, %.1f"  % (x,y)
            self._control.stage.goAbsolute(x,y)

    def offLasers(self):
        self._control.lasers.shutDown()
        self._control.illumination.enable561(False)
        print "Lasers off."

    def newSettings(self, settings):
        if type(settings) is str:
            if os.path.exists(settings):
                settingsDict = inLib.load_settings(settings)
                settingsDict['settings_filename'] = settings
                settings = settingsDict
        if type(settings) is dict:
            if not settings.has_key('settings_filename'):
                settings['settings_filename'] = 'unknown'
            self.savedSettings = [settings] + self.savedSettings
        return self.savedSettings

    def updateSettings(self, settings):
        if type(settings) is dict:
            if settings.has_key('devices'):
                self._control.camera.loadSettings(settings['devices']['camera']['settings'])
            elif settings.has_key('camera'):
                self._control.camera.loadSettings(settings['camera']['settings'])
            else:
                self._control.camera.loadSettings(settings)

    ###Talking to camera:

    def setFanMode(self, fanmode):
        self._control.camera.setFanMode(fanmode)

    def setExposureTime(self, exptime):
        self._control.camera.setExposureTime(exptime)

    def getExposureTimeFrameRate(self):
        self._control.camera._getAcquisitionTimings()
        fps = self._control.camera.getFrameRate()
        expTime = self._control.camera.getExposureTime()
        return expTime, fps

    def getROI(self):
        x0,x1,y0,y1 = self._control.camera.getROI()
        xres,yres = self._control.camera.getDimensions()
        self.dimensions = xres,yres
        return x0,y0,xres,yres

    def setROI(self, x0, y0, xres, yres, hbin, vbin):
        x1 = x0+xres-1
        y1 = x0+yres-1
        roi = [x0,x1,y0,y1]
        binning = [hbin, vbin]
        self._control.camera._setROIAndBinning(roi,binning)

    def setROI_commonsettings(self, settings):
        if settings == "128x128":
            self.setROI(193,193,128,128,1,1)
        if settings == "256x256":
            self.setROI(129,129,256,256,1,1)
        if settings == "512x512":
            self.setROI(1,1,512,512,1,1)
        if settings == "512x256":
            self.setROI(1,129,512,256,1,1)

    def openShutter(self):
        self._control.camera.openShutter()
        self.openShutterState=True

    def closeShutter(self):
        self._control.camera.closeShutter()
        self.openShutterState=False

    def getShutterState(self):
        return self.openShutterState

    def getGain(self):
        gain = self._control.camera.getEMCCDGain()
        return gain

    def setGain(self, gain):
        self._control.camera._setEMCCDGain(gain)

    ###    Acquiring data

    def getImageForPreview(self, modNum=None):
        if not self.filming:
            im = self._control.camera.getMostRecentImageNumpy()
            return im
        elif len(self.frames)==0:
            return None
        else:
            if self.frameToDisplay>=len(self.frames):
                self.frameToDisplay = 0
            if modNum is not None:
                if np.mod(self.zerothFrame-modNum[0],modNum[1])>0:
                    return None
            imBuffer = self.frames[self.frameToDisplay]
            im = np.fromstring(imBuffer,np.dtype('<u2'))
            return im.reshape(self._control.camera._settings['dimensions'])

    def setFixedLength(self, fixedLengthNum):
        self.fixedLengthNum = fixedLengthNum

    def useFixedLength(self, state):
        self.fixedLength = state

    def bufferStats(self):
        bufferSize = self._control.camera.getSizeOfCircBuffer()
        firstIm, lastIm = self._control.camera.getNumAvailableImages()
        return bufferSize, firstIm, lastIm

    def beginRecording(self, fileName):
        camSettings = self._control.camera.getSettings()
        self.daxfile = writer.DaxFile(fileName, camSettings)
        self.frames = []
        self.filming = 1
        self.recordedFrames = 0
        self.frameToDisplay = 0

    def recordToDAX(self):
        [frames, state] = self._control.camera.getImages16()
        done = False
        if state == "acquiring":
            if len(frames)>0:

                '''
                if len(self.frames)>0:
                    for i in range(len(frames)):
                        self.frames.extend([0])
                else:
                    self.frames = frames
                '''
                self.frames = frames

                if self.fixedLength:
                    toget = self.fixedLengthNum - self.recordedFrames
                    tosave = min(toget, len(frames))
                else:
                    tosave = len(frames)

                self.zerothFrame = self.recordedFrames

                if self.filming:
                    for i in range(tosave):
                        self.daxfile.saveFrames(frames[i],1,LEtoBE=True)
                        self.recordedFrames += 1
            
        else:
            print " run " + state

        if self.fixedLength:
            if self.recordedFrames >= self.fixedLengthNum:
                done = True
                self.filming = 0

        return self.frames, self.recordedFrames, done

    def endRecording(self):
        self.daxfile.closeFile([0,0,0],0)
        print " closing DAX file " + str(self.daxfile.filename)
        self.filming = 0

    def stopCamera(self):
        self._control.camera.stopCapture()

    def startCamera(self):
        self._control.camera.beginAcquisition()


            
            
        
        

