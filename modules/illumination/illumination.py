import inLib
import time
import numpy as np
import os, glob


class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing illumination..."
        inLib.Module.__init__(self, control, settings)
        print "illumination module initialized."

        self.initSettings = settings
        self.saved_settings = []
        #self.setSettings(settings)

        self.channelLED = 4
        self.channel561 = 3
        self.channelLED_AnalogOut = 1
        self.channel561_AnalogOut = 0

        #self.channelsA = [642, 488, 401, 561, 'LED']
        #self.channelsB = [0, 1, 2, 3, 4]

        self.laser_lines = np.array(self._control.lasers.laser_lines)
        
        w401 = np.where(abs(self.laser_lines - 401) < 10)
        w488 = np.where(abs(self.laser_lines - 488) < 10)
        w642 = np.where(abs(self.laser_lines - 642) < 10)

        self.laser_ports = [w401[0], w488[0], w642[0]]

        self._powerInMW = self._control.lasers._powerInMW

        #DAQ is self._control.DAQ
        #obis is self._control.lasers

        self.led_pow_func = lambda power: 5.0*power
        power_cal = np.array([2.53472938, -4.06563764, -0.15894481,
                              3.65015946, -2.4570782, 0.90482455,
                              0.02448929])
        self.sapphire_pow_func = lambda power: np.polyval(power_cal, power)

        self.shuttersRunning = False
        self.shuttersReady = False

    def newSettings(self, filename):
        settings_dict = inLib.load_settings(filename)
        if settings_dict.has_key('devices'):
            settings_dict = settings_dict['devices']['shutters']['settings']#['channel_settings']
        settings_dict['settings_filename'] = filename
        self.saved_settings = [settings_dict] + self.saved_settings
        return self.saved_settings

    def updateSettings(self, settings):
        self._control.shutters.stop()
        self.setShutters(settings['channels'],
                         settings['frames'],
                         settings['oversampling'],
                         settings['channel_settings'])
                         

    def enableLED(self, on):
        if on:
            self._control.shutters.digitalOutput(self.channelLED,1)
        else:
            self._control.shutters.digitalOutput(self.channelLED,0)

    def enable561(self, on):
        if on:
            self._control.shutters.digitalOutput(self.channel561,1)
        else:
            self._control.shutters.digitalOutput(self.channel561,0)

    def enableDiode(self, port, on):
        self._control.lasers.setLaserOnOff(port,on)

    def setDigitalMod(self, port, state):
        if state:
            self._control.lasers.setExtControl(port, "DIG")
        else:
            self._control.lasers.setInternalCW(port)

    def setLED(self, value):
        #Assuming LED on channel 1 Analog Out
        self._control.shutters.startAnalogOuts(self.channelLED_AnalogOut,value)

    def set561(self, value):
        #Assuming 561 on channel 0 Analog Out
        self._control.shutters.startAnalogOuts(self.channel561_AnalogOut,value)

    def setDiode(self, port, value):
        self._control.lasers.setPower(port, value)

    def setShutters(self, numChannels, frames, oversampling, channelSettings, freq=1000, noCounter=False):
        frametime = self._control.camera.getExposureTime(use_kinetic=True)
        waveforms = self._control.shutters.makeWaveform(numChannels,
                                                        frames,
                                                        oversampling,
                                                        channelSettings)
        self._control.shutters.readyShutters(frametime, frames, oversampling,
                                             numChannels, waveforms, freq=freq, noCounter=noCounter)
        self.shuttersReady = True
        return waveforms

    def ensureDigitalModCorrect(self, channelSettings):
        laserPorts = [self.laser_ports[2], self.laser_ports[1], self.laser_ports[0]]
        shouldBeDigMod = [False, False, False]
        for i in range(3):
            if (channelSettings[i]['ch_on'] > 0) or (channelSettings[i]['ch_off'] > 0):
                shouldBeDigMod[i] = True
                self.setDigitalMod(laserPorts[i], True)
        return shouldBeDigMod

    def shuttersUnready(self):
        self.shuttersReady = False

    def shuttersStop(self):
        self.shuttersRunning = False
        self._control.shutters.stop()
        self.resetDigitalOut()

    def startShutters(self, noCounter=False):
        if not self.shuttersRunning:
            self._control.shutters.start(noCounter=noCounter)
            self.shuttersRunning = True
        else:
            self.shuttersStop()
        return self.shuttersRunning

    def resetDigitalOut(self):
        for i in range(0,6):
            self._control.shutters.digitalOutput(i,0)
