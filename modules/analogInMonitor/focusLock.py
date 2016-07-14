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
        
        #self._control.shutters ==> controls NIDAQ


    def reInit(self):
        self.all_data = []
        self.stage_zs = []
        self.times = []
        self.means = []

        
    def getData(self, respond=False):
        data = self._control.shutters.getAI()
        self.data = np.frombuffer(data)
        data_mean = self.data.mean()
        if respond:
            current_z = self._control.piezo.getPosition(3)
            move_z = (data_mean - self.lock) * self.prop
            self.stage_zs.append(current_z)
            self._control.piezo.moveTo(3, current_z+move_z,
                                       waitForConvergence=False)
        self.resetMonitor()
        #self.all_data.append(self.data)
        self.means.append(data_mean)
        self.times.append(time.clock())
        return self.data, data_mean

    def getDataMeans(self):
        return self.means

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

    def startMonitor(self):
        self._control.shutters.configureAI()
        self._control.shutters.startAI()

    def stopMonitor(self):
        self._control.shutters.stopAI()

    def resetMonitor(self):
        self._control.shutters.stopAI()
        self._control.shutters.createAI()

    def saveData(self):
        #np.save("focusLockData.npy", self.all_data)
        np.save("focusLockDataMeans.npy", self.means)
        np.save("stage_zs.npy", self.stage_zs)
        np.save("times.npy", self.times)
        self.reInit()
