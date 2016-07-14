#!/usr/bin/python


import numpy as np
import time
import inLib


class Control(inLib.Module):

    def __init__(self, control, settings):
        print 'Initializing sliScan.'
        inLib.Module.__init__(self, control, settings)
        print 'sliScan initialized.'


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


    def calculateSLIScan(self, start, end, nSteps, filename=None):
        '''
        Calculates the modulation patterns for an  SLI Scan.

        :Parameters:
            *start*: float
                Scan start distance to the mirror in micrometers.
            *end*: float
                Scan end distance to the mirror in micrometers.
            *nSteps*: int
                The number of scan steps.
            *filename*: str
                The name of the file in which to store the scan patterns.

        :Returns:
            *scan_patterns*: numpy.array
                The scan patterns.
        '''
        n = self._control.slm.getNPixels()
        zs = np.linspace(start, end, nSteps)
        exp = self._control.slm.getSLIExperiment()
        dmf = self._control.slm.getSLIParams()['dmf']
        scan_patterns = np.zeros((nSteps,n,n))
        for i in xrange(nSteps):
            print 'sliscan: Calculating pattern {0}/{1}.'.format(i+1, nSteps)
            scan_patterns[i] = -np.angle(exp.get_sli_pupil_function(zs[i], dmf))
        if filename:
            np.save(filename, scan_patterns)
        return scan_patterns


    def scan(self, scan_file, nFrames=1, filename=None):
        '''
        Loads modulation patterns from a scan file iteratively.

        :Parameters:
            *scan_file*: str
                The file containing the modulation patterns.
            *filename*: str
                If not *None*, camera data for each pattern will be saved in this file
                in the.npy format.
            *nFrames*: int
                The number of camera frames to be taken for each pattern.
        '''
        dim = self._control.camera.getDimensions()
        scan = np.load(scan_file, mmap_mode='r')
        data = np.zeros((scan.shape[0],) + dim)
        slicesFrames = np.zeros((nFrames,)+dim)
        frame_length = 1.0/self._control.camera.getFrameRate()
        self._control.slm.setActive(True)
        self._control.slm.setOthersActive(True)
        index = self._control.slm.addOther(np.zeros_like(scan[0]))
        for i in xrange(scan.shape[0]):
            print 'sliscan: Scanning step {0}/{1}.'.format(i+1, scan.shape[0])
            self._control.slm.setOtherData(index, scan[i])
            time.sleep(2*frame_length)
            for j in xrange(nFrames):           
                slicesFrames[j] = self._control.camera.getMostRecentImageNumpy()
                time.sleep(frame_length)
            data[i] = np.mean(slicesFrames, axis=0)
        self._control.slm.deleteOther(index)
        if filename:
            np.save(filename, data)
        return data
