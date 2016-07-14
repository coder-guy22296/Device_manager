import inLib
import time
import numpy as np
from scipy.fftpack import fft2,ifft2,fftshift
import os, glob
import serial


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

        #self._control.xystage ==> asi

        #self._control.syringePumps 

        self.numSavedLocs = 10

        self.savedLocations = []
        self.savedSettings = [self.initSettings]

        #For arduino (can shutter [slowly] laser)
        self.baudrate=57600
        self.port = "COM18"

        try:
            self.ser = serial.Serial(self.port, self.baudrate)
            serin = self.ser.read()
            print "Arduino connection: ", serin
        except:
            print "Unable to make Arduino connection."
            self.ser = None


    def reInit(self):
        self.all_data = []
        self.stage_zs = []
        self.times = []
        self.means = []

    def stopPump(self, num):
        self._control.syringePumps.stop(num)

    def runPumpForward(self,num):
        self._control.syringePumps.runForward(num)

    def setPumpForwardRate(self,num, flowrate):
        #Flowrate in ulh
        self._control.syringePumps.setForwardRate(num,flowrate,0)

    def writePositionToArduino(self, pos):
        if self.ser is not None:
            self.ser.write("%i,%i,%i,%i" % (1,pos,0,1))

    def getXYZPosition(self):
        x,y = self._control.xystage.whereStage()
        z = None
        return x,y,z

    def getStageStatus(self):
        resp = self._control.xystage.getStatus()
        #print "Length of status: ", len(resp)
        #print "Modules response for stage status: ", resp[0]=='N'
        #return resp[0]=='N'
        #print "Stage resp: ", resp
        if resp[0]=='M' or resp[1]=='M':
            #print "Stage stat False"
            return False
        else:
            #print "Stage stat True"
            return True

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
        #self.savedLocations = []
        self._control.xystage.clearLocations()

    def addLocation(self):
        locs = self._control.xystage.addLocation()
        return locs

    def savePositions(self, filename):
        self._control.xystage.saveLocations(filename)

    def loadPositions(self, filename):
        locs = self._control.xystage.loadLocations(filename)
        return locs

    def goToLocations(self,num):
        self._control.xystage.goToLocations(num)

    def setWaitTime(self, wtime):
        self._control.xystage.setWaitTime(wtime)

    def setSpeed(self, axis, speed):
        self._control.xystage.setSpeed(axis, speed)
        
    def returnWaitTime(self):
        wtime = self._control.xystage.returnWaitTime()
        return wtime
        
    '''
    def goToLocation(self, num):
        if num<=len(self.savedLocations):
            x,y = self.savedLocations[num-1]
            print "Moving to %.1f, %.1f"  % (x,y)
            self._control.xystage.goAbsoluteXY(x,y)
    '''

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
