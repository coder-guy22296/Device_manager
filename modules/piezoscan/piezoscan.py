#!/usr/bin/python


import numpy as np
import time
import inLib


class Control(inLib.Module):

    def __init__(self, control, settings):
        print 'Initializing Piezoscan.'
        inLib.Module.__init__(self, control, settings)
        print 'Piezoscan initialized.'
        self.active = False

        if settings['ThorlabsMotor'] == True:
            self.useThorlabs = True
        else:
            self.useThorlabs = False

    '''
    @property
    def active(self):
        return self._active

    @active.setter
    def active(self, value):
        print 1
        if type(value) == bool:
            print 2
            self._active = value
            print 3
        else:
            raise TypeError
    '''

    def calcScanParams(self, start, end, nSteps):
        '''
        Calculates the step size of a scan.
        
        :Parameters:
            *start*: float
                Start of scan in micrometers, relative to current position.
            *end*: float
                End of scan in micrometers, relative to current position.
            *nSteps*: int
                The number of steps.

        :Returns:
            *stepSize*: float
                The step size in nanometers.
        '''
        stepSize = abs(1000.0*(end-start)/(nSteps-1.0))
        return stepSize

    def scan(self, start, end, nSteps, nFrames, filename=None):
        '''
        Performs a scan
        '''
        if self.useThorlabs:
            if start<end:
                up=True
            else:
                up=False
            data = self.scan_thorlabs(nSteps, nFrames, up=up, filename=filename)
        else:
            data = self.scan_piezo(start, end, nSteps, nFrames, filename=filename)
        return data


    def scan_thorlabs(self, nSteps, nFrames, up=True, filename=None):
        self.active = True
        dim = self._control.camera.getDimensions()
        data = np.zeros((nSteps,) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        for i in xrange(nSteps):
            if self.active:
                if up:
                    self._control.servo.jogUp()
                else:
                    self._control.servo.jogDown()
                time.sleep(4*frame_length)
                for j in xrange(nFrames):
                    im = self._control.camera.getMostRecentImageNumpy()
                    if im is None:
                        time.sleep(frame_length)
                        im = self._control.camera.getMostRecentImageNumpy()
                    slicesFrames[j] = im
                    time.sleep(frame_length)
                data[i] = np.mean(slicesFrames, axis=0)
            else:
                break
        if self.active and filename:
            print 'piezoscan: Saving scan to', filename
            np.save(filename, data)
            self.active = False
        return data

    def scan_piezo(self, start, end, nSteps, nFrames, filename=None):
        '''
        Performs a piezo scan.

        :Parameters:
            *start*: float
                The starting point relative to the current position, in micrometers.
            *end*: float
                The end point, relative to the current position, in micrometers.
            *nSteps*: int
                The number of steps.
            *nFrames*: int
                The number of frames to be averaged in each step.
            *filename*: str
                If not *None*, the image data will be saved in a .npy file with this name.
        '''
        print 'piezoscan: Scanning with params:', start, end, nSteps
        self.active = True
        dz = abs((end-start)/(nSteps-1.0))
        orig_z = self._control.piezo.getPosition(3)
        start += orig_z
        end += orig_z
        zs = np.linspace(start, end, nSteps)
        dim = self._control.camera.getDimensions()
        data = np.zeros((nSteps,) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        for i in xrange(nSteps):
            if self.active:
                self._control.piezo.moveTo(3, zs[i])
                time.sleep(2*frame_length)
                for j in xrange(nFrames):           
                    slicesFrames[j] = self._control.camera.getMostRecentImageNumpy()
                    time.sleep(frame_length)
                data[i] = np.mean(slicesFrames, axis=0)
            else:
                break
        self._control.piezo.moveTo(3, orig_z)
        if self.active and filename:
            print 'piezoscan: Saving scan to', filename
            np.save(filename, data)
            self.active = False
        return data

    def stop(self):
        self.active = False
