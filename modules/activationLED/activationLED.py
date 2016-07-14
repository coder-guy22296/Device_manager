import inLib
import time
import numpy as np
import os, glob


class Control(inLib.Module):

    def __init__(self, control, settings):
        print "Initializing activationLED..."
        inLib.Module.__init__(self, control, settings)
        print "activationLED module initialized."

        self.initSettings = settings
        self.saved_settings = []
        
        self.channelLED = 4
        self.channelLED_AnalogOut = 1

        self.led_pow_func = lambda power: 5.0*power

        self.running = False
        self.ready = False

    def updateSettings(self, settings):
        return 0

    def setLED(self, value):
        #Assuming LED on channel 1 Analog Out
        self._control.shutters.startAnalogOuts(self.channelLED_AnalogOut,value)
                         

    def enableLED(self, on):
        if on:
            self._control.shutters.digitalOutput(self.channelLED,1)
        else:
            self._control.shutters.digitalOutput(self.channelLED,0)

    def resetDigitalOut(self):
        for i in range(0,6):
            self._control.shutters.digitalOutput(i,0)
