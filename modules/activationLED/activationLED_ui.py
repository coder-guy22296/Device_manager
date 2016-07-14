#!/usr/bin/python

from PyQt4 import QtCore,QtGui,Qwt5
from functools import partial
import inLib
import sys, time
import numpy as np

from Utilities import QExtensions as qext
 
class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        design_path = 'modules.activationLED.activationLED_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        self._ui.OnLED.stateChanged.connect(partial(self._onThroughDAQ,1))
        self._ui.SliderLED.valueChanged.connect(partial(self._powerThroughDAQ,1))

        self._ui.pushButton_start.setEnabled(True)
        self._ui.pushButton_start.clicked.connect(self._start)

        self._activator = None
        self._repeater = QtCore.QTimer()
        self._repeater.timeout.connect(self._repeatCycles)
        
        self.settings = []

        self.cycles = 0
        self.numCycles = 100
        


    def _newSettings(self, filename):
        self.settings = self._control.newSettings(filename)
        print self.settings
        if len(self.settings) > self._max_saved:
            self.settings.pop()
        for i,p in enumerate(self.settings):
            filename = p["settings_filename"]
            radiobutton = getattr(self._ui, "settings"+str(i+2)+"_radioButton")
            radiobutton.setText(filename.split('/')[-1][:-5])
            radiobutton.show()
        self._ui.settings2_radioButton.click()

    def _repeatCycles(self):
        if self._activator is not None:
            self._activator.start()
            print "Started activator thread..."
            self.cycles += 1
        if self.cycles == self.numCycles:
            self._stop()

    def _start(self):
        if self._ui.pushButton_start.text() == 'Start':
            self.cycles = 0
            onTime = float(self._ui.lineEdit_onTime.text())
            periodTime = float(self._ui.lineEdit_periodTime.text())
            self.numCycles = int(self._ui.lineEdit_cycles.text())
            self._activator = Activate(self._control, onTime)
            self._ui.pushButton_start.setText('Stop')
            self._repeater.start(periodTime*1000)
            #self._activator.start()
        else:
            self._stop()

    def _stop(self):
        self._repeater.stop()
        self._ui.pushButton_start.setText('Start')
        self._activator = None
            
    def _onThroughDAQ(self, channel):
        if channel:
            state = self._ui.OnLED.isChecked()
            self._control.enableLED(state)

    def _powerThroughDAQ(self, channel):
        if channel:
            power = self._ui.SliderLED.value()
            voltage = self._control.led_pow_func(power/100.0)
            self._control.setLED(voltage)


class Activate(QtCore.QThread):
    def __init__(self, control, onTime):
        QtCore.QThread.__init__(self)
        self.control = control
        self.onTime = onTime

    def run(self):
        self.control.enableLED(True)
        print "LED is on..."
        time.sleep(self.onTime)
        self.control.enableLED(False)
        print "LED is off..."

    def stop(self):
        self.control.enableLED(False)
        return 0
        


