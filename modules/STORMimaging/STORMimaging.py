import inLib
import time
import numpy as np
from scipy.fftpack import fft2,ifft2,fftshift
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

        self.numSavedLocs = 8

        self.savedLocations = []
        self.savedSettings = [self.initSettings]


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
        if len(self.savedLocations)<self.numSavedLocs:
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

    def setExposureTime(self, exptime):
        self._control.camera.setExposureTime(exptime)

    def getExposureTimeFrameRate(self):
        #self._control.camera._getAcquisitionTimings()
        fps = self._control.camera.getFrameRate()
        expTime = self._control.camera.getExposureTime()
        return expTime, fps

    def getROI(self):
        x0,x1,y0,y1 = self._control.camera.getROI()
        xres,yres = self._control.camera.getDimensions()
        self.dimensions = xres,yres
        return x0,y0,xres,yres
