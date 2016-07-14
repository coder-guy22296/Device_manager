#!/usr/bin/python


import inLib
from PyQt4.QtGui import QImage
import numpy as np
import libs.imagewriters as writer
import time
import ctypes
import Queue
import threading
#from modules.adaptiveOptics import signalForAO

prop_ids = {"bits per channel" : 4325680,
            "trigger source" : 1048848,
            "trigger mode" : 1049104,
            "trigger active" : 1048864,
            "trigger polarity" : 1049120,
            "trigger connector" : 1049136,
            "trigger times" : 1049152,
            "trigger delay" : 1049184,
            "exposure time" : 2031888,
            "defect correct mode" : 4653072,
            "binning" : 4198672,
            "subarray hpos" : 4202768,
            "subarray hsize" : 4202784,
            "subarray vpos" : 4202800,
            "subarray vsize" : 4202816,
            "subarray mode" : 4202832,
            "timing readout time" : 4206608,
            "timing cyclic trigger period" : 4206624,
            "timing min trigger blanking" : 4206640,
            "timing min trigger interval" : 4206672,
            "timing exposure" : 4206688,
            "timing invalid exposure period" : 4206704,
            "internal frame rate" : 4208656,
            "internal frame interval" : 4208672,
            "image width" : 4325904,
            "image height" : 4325920,
            "image rowbytes" : 4326192,
            "image framebytes" : 4325952,
            "image top offset bytes" : 4325968,
            "image pixel format" : 4326000,
            "buffer rowbytes" : 4326192,
            "buffer framebytes" : 4326208,
            "buffer top offset bytes" : 4326224,
            "buffer pixel type" : 4326240,
            "number of output trigger connector" : 1835024,
            "output trigger polarity[0]" : 1835296,
            "output trigger active[0]" : 1835312,
            "output trigger delay[0]" : 1835328,
            "output trigger period[0]" : 1835344,
            "output trigger kind[0]" : 1835360}

cam_info_strings = {"Model" : 0x04000104,
                    "CameraID" : 0x04000102,
                    "Camera Version" : 0x04000105,
                    "Driver Version" : 0x04000106,
                    "Module Version" : 0x04000107,
                    "DCAM-API Version" : 0x04000108}

'''Trigger modes
/*** --- trigger mode --- ***/
enum {
	DCAM_TRIGMODE_INTERNAL				= 0x0001,
	DCAM_TRIGMODE_EDGE					= 0x0002,
	DCAM_TRIGMODE_LEVEL					= 0x0004,
	DCAM_TRIGMODE_MULTISHOT_SENSITIVE	= 0x0008,
	DCAM_TRIGMODE_CYCLE_DELAY			= 0x0010,
	DCAM_TRIGMODE_SOFTWARE				= 0x0020,
	DCAM_TRIGMODE_FASTREPETITION		= 0x0040,
	DCAM_TRIGMODE_TDI					= 0x0080,
	DCAM_TRIGMODE_TDIINTERNAL			= 0x0100,
	DCAM_TRIGMODE_START					= 0x0200,
	DCAM_TRIGMODE_SYNCREADOUT			= 0x0400
	, _DCAM_TRIGMODE_MULTIGATE			= 0x0800	/*[[ internal ]]*/
	, _DCAM_TRIGMODE_MULTIGATE_LEVEL	= 0x1000	/*[[ internal ]]*/

'''
            
            

class Control(inLib.Device):
    '''
    
    '''

    def __init__(self, settings):
        inLib.Device.__init__(self, 'hamamatsu.orcaflash.camera_api', settings)

        self._api.openCamera()

        #Camera properties
        self._props = {}

        #Loaded yaml settings
        self._initSettings = settings.copy()
        self.loadedSettings = settings
        self.saved_settings = []
        self.newSettings('no_file', settings_dict=settings)
        self.updateSettings(settings)

        self._previewFrames = 10
        self._recordFrames = 100

        self.daxfile = 0

        self.recording=False

        self.fixedLength = False
        self.fixedLengthNum = 100

        self.latestFrame = -1
        self.j=0

        self.frame_bytes = 0
        self.frame_x = 0
        self.frame_y = 0

        #Setting output triggering
        self._setProp(1835360, 3) #Programmable mode
        self._setProp(1835344, 0.001) #1ms period
        self._setProp(1835296, 2) #positive polarity

    def printProperties(self):
        nl = '\n'
        for key in prop_ids:
            value = self._api.getPropertyValue(prop_ids[key])
            print key + str(": ") + str(value) + nl

    def setOutputTriggerPeriod(self, period):
        self._setProp(1835344, period)
        print "New trigger period: ", self._getProp(1835344)

    def setToExposureMode(self, polarity=1):
        self._setProp(1835360, 2)
        self._setProp(1835296, polarity)

    def setToProgrammableMode(self, polarity=2, period=0.001):
        self._setProp(1835360, 3) #Programmable mode
        self._setProp(1835344, period) #1ms period
        self._setProp(1835296, polarity) #positive polarity
            

    def newSettings(self, filename, settings_dict=None):
        if settings_dict is None:
            settings_dict = inLib.load_settings(filename)
            if settings_dict.has_key('devices'):
                settings_dict = settings_dict['devices']['camera']['settings']
            print settings_dict
            settings_dict['settings_filename'] = filename
        self.saved_settings = [settings_dict] + self.saved_settings
        return self.saved_settings

    def updateSettings(self, settings_dict):
        try:
            roi = settings_dict['roi']
            x0 = roi[0] 
            y0 = roi[2] 
            xsize = roi[1] - x0
            ysize = roi[3] - y0
            print "Setting ROI to: %i,%i. And size: %i,%i." % (x0,y0,xsize,ysize)
            self._setROI(x0,y0,xsize,ysize)
            time.sleep(0.05)
            self._setROI(x0,y0,xsize,ysize)
        except:
            print "Unable to set new ROI"
        try:
            exposure_time = settings_dict['exposure_time']
            print "Setting exposure time to: %.2f ms" % (exposure_time*1000)
            self.setExposureTime(exposure_time, inMillisec=False)
        except:
            print "Unable to set new timings"

    def loadSettings(self, settings):
        if type(settings) is str:
            if os.path.exists(settings):
                settingsDict = inLib.load_settings(settings)
                settings = settingsDict
        if type(settings) is dict:
            try:
                roi = settings['roi']
                x0 = roi[0] 
                y0 = roi[2] 
                xsize = roi[1] - x0
                ysize = roi[3] - y0
                print "Setting ROI to: %i,%i. And size: %i,%i." % (x0,y0,xsize,ysize)
                self._setROI(x0,y0,xsize,ysize)
                time.sleep(0.05)
                self._setROI(x0,y0,xsize,ysize)
            except:
                print "Unable to set new ROI"
            try:
                exposure_time = settings['exposure_time']
                print "Setting exposure time to: %.2f ms" % (exposure_time*1000)
                self.setExposureTime(exposure_time, inMillisec=False)
            except:
                print "Unable to set new timings"

        

    def _getStatus(self):
        s = self._api.getStatus()
        if s==0:
            status = "Error"
        elif s==1:
            status = "Busy"
        elif s==2:
            status = "Ready"
        elif s==3:
            status = "Stable"
        elif s==4:
            status = "Unstable"
        return status

    def _set2x2Binning(self):
        self._api.setPropertyValue(4198672, 2)

    def _setNoBinning(self):
        self._api.setPropertyValue(4198672, 1)

    def setDefectCorrect(self, defectMode):
        if defectMode==1 or defectMode==2:
            self._api.setPropertyValue(4653072,defectMode)
         
    def _setROI(self, x0, y0, xsize, ysize):
        self._api.setPropertyValue(4202832,2)   #Allows use of ROI
        time.sleep(0.05)
        self._api.setPropertyValue(4202768, x0)
        time.sleep(0.05)
        self._api.setPropertyValue(4202800, y0)
        time.sleep(0.05)
        self._api.setPropertyValue(4202784, xsize)
        time.sleep(0.05)
        self._api.setPropertyValue(4202816, ysize)
        time.sleep(0.05)
        self._props['dimensions'] = self._api.getDataSize()
        x0,y0 = self._getX0Y0()

    def setROI_commonsettings(self, setting):
        print "Common settings: ", setting
        if setting == "128x128":
            #self._setROI(960,960,128,128)
            self._setROI(1088,960,128,128)
        elif setting == "256x256":
            #self._setROI(768,896,256,256)
            self._setROI(1024,896,256,256)
        elif setting == "512x256":
            #self._setROI(768,896,512,256)
            self._setROI(768,896,512,256)
        elif setting == "512x512":
            #self._setROI(768,768,512,512)
            self._setROI(768,768,512,512)

    def _setProp(self, propID, value):
        self._api.setPropertyValue(propID, value)

    def _getProp(self, propID):
        val = self._api.getPropertyValue(propID)
        return val

    def _getX0Y0(self):
        x0 = self._api.getPropertyValue(4202768)
        y0 = self._api.getPropertyValue(4202800)
        self._props['x_start'] = x0
        self._props['y_start'] = y0
        return x0,y0

    def _getResolution(self):
        xres, yres = self._api.getDataSize()
        self._props['dimensions'] = xres,yres
        return xres, yres

    def getDimensions(self):
        xres, yres = self._api.getDataSize()
        self._props['dimensions'] = xres,yres
        return xres, yres

    def getROI(self):
        x0,y0 = self._getX0Y0()
        xres,yres = self._getResolution()
        return x0,x0+xres,y0,y0+yres

    #def _getExposureTime(self):
    #    expTime = self._api.getExposureTime()
    #    self._props['exposure_time'] = expTime
    #    return expTime

    def setTriggerMode(self, value):
        self.stopCapture()
        self._api.setTriggerMode(value)
        self.beginPreview()
        cyclicTriggerPeriod = self._getProp(prop_ids["timing cyclic trigger period"])
        print "Cyclic trigger period (sec): %.3f\n" % cyclicTriggerPeriod

    def getExposureTime(self, use_kinetic=False):
        expTime = self._api.getExposureTime()
        self._props['exposure_time'] = expTime
        return expTime

    def getFrameRate(self):
        internalFrameRate = self._api.getInternalFrameRate()
        self._props['internal_frame_rate'] = internalFrameRate
        return internalFrameRate

    #def _setExposureTime(self,exposure, inMillisec=False):
    #    if inMillisec:
    #        exposure *= 1e-3
    #    self._api.setExposureTime(exposure)


    def setExposureTime(self,exposure, inMillisec=False):
        if inMillisec:
            exposure *= 1e-3
        self._api.setExposureTime(exposure)

    def _continuousCapture(self, frames):
        self._api.preCapture(1)
        self._api.allocFrame(frames)
        self._api.capture()

    def saveFrameToNumpy(self, filename):
        xres,yres = self._props['dimensions']
        newestFrame = self._api.getTransferInfo()[0]
        if newestFrame>=0:
            im = self._api.lockData(newestFrame, xres*yres)
            np.save(filename, im.reshape(xres,yres))

    def getMostRecentImageNumpy(self):
        xres,yres = self._props['dimensions']
        newestFrame = self._api.getTransferInfo()[0]
        if newestFrame>=0:
            im = self._api.lockData(newestFrame, xres*yres)
            return im.reshape(xres,yres)
        else:
            return None

    def beginPreview(self):
        print "Starting preview mode..."
        xres,yres = self._getResolution()
        self._props['dimensions'] = xres,yres
        self._continuousCapture(self._previewFrames)

    def getImageForPreview(self):
        xres,yres = self._props['dimensions']
        newestFrame = self._api.getTransferInfo()[0]
        if newestFrame>=0:
            prev_im = self._api.lockData(newestFrame, xres*yres)
        else:
            return None
        self._api.unlockData()
        if prev_im is not None:
            return prev_im.reshape(xres,yres)
        else:
            return None

    def _setRecordMode(self, mode, length):
        self.fixedLength = bool(mode)
        self.fixedLengthNum = length

    def getLatest(self):
        xres,yres = self._props['dimensions']
        previousFrame = self.latestFrame
        totalFrames=0
        newestFrame, totalFrames = self._api.getTransferInfo()
        while totalFrames==0:
            #time.sleep(50)
            newestFrame, totalFrames = self._api.getTransferInfo()
        #print "newest, total: ", (newestFrame, totalFrames)
        if previousFrame==-1:
            numFrames = newestFrame
        else:
            numFrames = newestFrame - previousFrame
        if numFrames < 0:
            numFrames = newestFrame + (self._recordFrames - previousFrame)
        self.ims = np.zeros((numFrames, xres*yres),dtype=np.uint16)
        #ims = []
        for self.j in range(0,numFrames):
            frameToGet = previousFrame + 1 + self.j
            if frameToGet >= self._recordFrames:
                #print "looped..."
                frameToGet = frameToGet - self._recordFrames
            #print "Frame to Get: ", frameToGet
            #print "Size of self.ims: ", np.size(self.ims)
            self.ims[self.j] = self._api.lockData(frameToGet, xres*yres)
            #ims.append(self._api.lockData(frameToGet,xres*yes).tostring())
        self._api.unlockData()
        self.latestFrame = newestFrame
        return totalFrames

    def imageForPreviewWhileRecording(self):
        xres,yres = self._props['dimensions']
        if self.j>0 and len(self.ims)>self.j:
            return self.ims[self.j].reshape(xres,yres)
        
    def stopCapture(self):
        self._api.idle()
        self._api.freeFrame()      

    def beginRecording(self, fileName):
        self.stopCapture()
        self.daxfile = writer.DaxFile(fileName, self._props)
        self.latestFrame = -1
        self.recordedFrames = 0
        self.recording = True
        print "Starting to record..."
        self._continuousCapture(self._recordFrames)
        '''
        while totalFrames<numFrames and self.recording:
            time.sleep(5*self._props['exposure_time'])
            ims, latestFrame, totalFrames = self.getLatest(latestFrame)
            for i in range(0,len(ims)):
                self.daxfile.saveFrames(ims[i], 1)
        '''

    def beginSnaps(self, fileName):
        self.stopCapture()
        self.daxfile = writer.DaxFile(fileName, self._props)
        self.latestFrame = -1
        self.recordedFrames = 0
        self.recording = True
        print "Starting to record..."
        self._continuousCapture(self._recordFrames)

    def takeSnap(self):
        self.getLatest()
        #print "For snaps: len(ims), j: ", (len(self.ims),self.j)
        if self.j>0 and len(self.ims)>self.j:
            im = self.ims[self.j]
            temp = im.astype(np.dtype('>H')).copy()
            self.daxfile.saveFrames(temp.tostring(), 1)
            self.recordedFrames+=1
        return self.recordedFrames

    def recordToDAX(self):
        totalFrames = self.getLatest()
        num_ims = self.ims.shape[0]
        done = False
        if self.fixedLength:
            toget = self.fixedLengthNum - self.recordedFrames
            tosave = min(toget, num_ims)
        else:
            tosave = num_ims
        for i in range(0,tosave):
            temp = self.ims[i].astype(np.dtype('>H')).copy()
            self.daxfile.saveFrames(temp.tostring(), 1)
            self.recordedFrames += 1
            if self.fixedLength:
                if self.recordedFrames == self.fixedLengthNum:
                    done=True
        return totalFrames, self.recordedFrames, done
        

    def endRecording(self):
        #self._api.wait(4, 0x80000000)
        self._api.idle()
        time.sleep(1)
        self.daxfile.closeFile([0,0,0], 0)
        print "Closing DAX file " + str(self.daxfile.filename)
        self.latestFrame = -1

    def shutDown(self):
        print 'Shutting down ORCA Flash4. camera...'
        self.stopCapture()
        self._api.closeCamera()

