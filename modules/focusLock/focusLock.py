import inLib
import time
import numpy as np
from scipy.fftpack import fft2,ifft2,fftshift
import os, glob


class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing focusLock..."
        inLib.Module.__init__(self, control, settings)
        print "focusLock initialized."

        self.initSettings = settings
        
        self.setSettings(settings)
        self.reInit()

        self.prop = 0.1
        self.lock = -1.0

        self.offset = 0.000

        self.oneAIChannel = True

        self.intensityThreshold = 0
        
        #self.piezo = self._control.piezo
        #self._control.shutters ==> controls NIDAQ


    def reInit(self):
        self.all_data = []
        self.stage_zs = []
        self.times = []
        self.means = []
        self.means2 = []

    def setPropFeedback(self, prop):
        self.prop = prop

    def setOffset(self, offsetValue):
        self.offset = offsetValue

    def setIntensityThreshold(self, thresh):
        self.intensityThreshold = thresh

    def findOffset(self):
        current_z = 0
        numIterations = 3
        for i in range(0,numIterations):
            self._control.piezo.moveTo(3,50.0)
            time.sleep(0.01)
            current_z += self._control.piezo.getPosition(3)
        offset = (current_z / numIterations) - 50.0
        return offset

    def setLock(self, lock):
        self.lock = lock
        
    def getData(self, respond=False):
        enoughIntensity = True
        data = self._control.shutters.getAI()
        self.data = np.frombuffer(data)
        if not self.oneAIChannel:
            dataLen = len(self.data)/2
            self.data2 = self.data[dataLen:]
            self.data = self.data[:dataLen]
        data_mean = self.data.mean()
        respondTo = data_mean
        if not self.oneAIChannel:
            data_mean2 = self.data2.mean()
            respondTo = data_mean / (2.0*data_mean2)
            if data_mean2 < self.intensityThreshold:
                enoughIntensity = False
        if respond and enoughIntensity:
            current_z = self._control.piezo.getPosition(3) - self.offset
            move_z = (respondTo - self.lock) * self.prop
            self.stage_zs.append(current_z)
            if abs(move_z)<10:
                self._control.piezo.moveTo(3, current_z+move_z,
                                           waitForConvergence=False)
        self.resetMonitor()
        #self.all_data.append(self.data)
        self.means.append(data_mean)
        if not self.oneAIChannel:
            self.means2.append(data_mean2)
        self.times.append(time.clock())
        if self.oneAIChannel:
            return self.data, data_mean
        else:
            return self.data, data_mean, data_mean2

    def getDataMeans(self):
        if self.oneAIChannel:
            return self.means
        else:
            return self.means, self.means2

    def setCurrentLock(self):
        if len(self.means)>2:
            self.lock = self.means[-1]
            return self.lock
        else:
            return False

    def setSettings(self, settings):
        self.samples = settings['samples']
        self.sampleRate = settings['sampleRate']
        self._control.shutters.setAIParams(self.samples, self.sampleRate)

    def calibrate(self, distance):
        self._control.shutters.setAIParams(500, 10000)
        self.startMonitor()
        time.sleep(1.0)
        current_z = self._control.piezo.getPosition(3)
        cal_distances = np.linspace(current_z-distance, current_z+distance,100)
        cal_values = np.zeros_like(cal_distances)
        i=0
        data, initial_value = self.getData(respond=False)
        self.startMonitor()
        for position in cal_distances:
            print i
            self._control.piezo.moveTo(3,position,waitForConvergence=False)
            time.sleep(0.05)
            data, cal_values[i] = self.getData(respond=False)
            i+=1
            self.startMonitor()
        self._control.piezo.moveTo(3,current_z,waitForConvergence=False)
        self.resetMonitor()
        fit = np.polyfit(cal_values[1:-1] - initial_value, cal_distances[1:-1], 1)
        fitvals = np.polyval(fit, cal_distances[1:-1])
        errors = abs(cal_values[1:-1] - initial_value - fitvals)/fitvals
        return cal_distances, cal_values, fit[0], errors

    def startMonitor(self):
        if not self.oneAIChannel:
            self._control.shutters.addAIChannel()
        self._control.shutters.configureAI()
        self._control.shutters.startAI()

    def stopMonitor(self):
        self._control.shutters.stopAI()

    def resetMonitor(self):
        self._control.shutters.stopAI()
        self._control.shutters.createAI()

    def useTwoChannels(self, state):
        self.oneAIChannel = not state

    def saveData(self):
        #np.save("focusLockData.npy", self.all_data)
        np.save("focusLockDataMeansCh1.npy", self.means)
        np.save("focusLockDataMeansCh2.npy", self.means2)
        np.save("stage_zs.npy", self.stage_zs)
        np.save("times.npy", self.times)
        self.reInit()
