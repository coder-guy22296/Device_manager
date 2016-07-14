import inLib
import time
import numpy as np
import os, glob


class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing testing..."
        inLib.Module.__init__(self, control, settings)
        print "illumination module initialized."

        self.initSettings = settings

    '''

    def enableLED(self, on):
        if on:
            self._control.DAQ.digitalOutput(self.channelLED,1)
        else:
            self._control.DAQ.digitalOutput(self.channelLED,0)

    def enable561(self, on):
        if on:
            self._control.DAQ.digitalOutput(self.channel561,1)
        else:
            self._control.DAQ.digitalOutput(self.channel561,0)

    def enableDiode(self, port, on):
        self._control.obis.setLaserOnOff(port,on)

    def setLED(self, value):
        #Assuming LED on channel 1 Analog Out
        self._control.DAQ.startAnalogOuts(1,value)

    def set561(self, value):
        #Assuming 561 on channel 0 Analog Out
        self._control.DAQ.startAnalogOuts(0,value)

    def setDiode(self, port, value):
        self._control.obis.setPower(port, value)
        '''
