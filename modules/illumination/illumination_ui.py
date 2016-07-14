#!/usr/bin/python

from PyQt4 import QtCore,QtGui,Qwt5
from functools import partial
import inLib
import sys
import numpy as np

from Utilities import QExtensions as qext
 
class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        design_path = 'modules.illumination.illumination_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        laser_lines = self._control.laser_lines
        laser_ports = self._control.laser_ports

        self.channelsA = [642, 488, 401, 561, 'LED']
        self.channelsB = [0, 1, 2, 3, 4]
        
        self.sliders = [self._ui.Slider401, self._ui.Slider488, self._ui.Slider642]
        self.onoff = [self._ui.On401, self._ui.On488, self._ui.On642]
        self.groupBoxes = [self._ui.OBIS401, self._ui.OBIS488, self._ui.OBIS642]
        self.powerTxts = [self._ui.Pow401, self._ui.Pow488, self._ui.Pow647]
        color = ["(135,0,205,80)", "(0,240,255,80)", "(255,10,0,80)"]
        self.digMod = [self._ui.checkBox_digMod401, self._ui.checkBox_digMod488, self._ui.checkBox_digMod642]

        for i in range(0,len(laser_ports)):
            if len(laser_ports[i])==1:
                print (i, laser_ports[i][0])
                self.sliders[i].valueChanged.connect(partial(self._updatePower, i, laser_ports[i][0]))
                self.onoff[i].stateChanged.connect(partial(self._turnOnOff, i, laser_ports[i][0]))
                self.digMod[i].stateChanged.connect(partial(self._setDigitalMod, i, laser_ports[i][0]))
                self.groupBoxes[i].setStyleSheet("QGroupBox {background-color: rgba"+color[i]+"}")
            else:
                self.sliders[i].setDisabled(True)
                self.onoff[i].setDisabled(True)
                self.groupBoxes[i].hide()

        self._ui.Sapphire561.setStyleSheet("QGroupBox {background-color: rgba(218,255,0,80)}")
        self._ui.LED.setStyleSheet("QGroupBox {background-color: rgba(255,255,255,80)}")

        self._ui.On561.stateChanged.connect(partial(self._onThroughDAQ,0))
        self._ui.OnLED.stateChanged.connect(partial(self._onThroughDAQ,1))
        self._ui.Slider561.valueChanged.connect(partial(self._powerThroughDAQ,0))
        self._ui.SliderLED.valueChanged.connect(partial(self._powerThroughDAQ,1))

        self._ui.checkBox_enableShuttering.stateChanged.connect(self._enableShuttering)
        self._ui.lineEdit_freq.setText("1")
        self._ui.lineEdit_shutterFrames.setText("10")
        self._ui.lineEdit_oversampling.setText("100")

        self.useShuttering = False
        self._ui.pushButton_startShutters.setEnabled(False)
        self._ui.pushButton_startShutters.clicked.connect(self._startShutters)
        self._ui.pushButton_readyShutters.setEnabled(False)
        self._ui.pushButton_readyShutters.clicked.connect(self._getShutterParams)

        self._ui.groupBox_shuttering401.hide()
        self._ui.groupBox_shuttering488.hide()
        self._ui.groupBox_shuttering561.hide()
        self._ui.groupBox_shuttering642.hide()
        self._ui.groupBox_shutteringLED.hide()

        self._ui.pushButton_loadFile.clicked.connect(self.loadShutterFile)

        self._ui.checkBox_561.hide()
        self._ui.checkBox_LED.hide()

        #handling file drops:
        
        #self._ui.widget_shuttering.dragEnterEvent = self.dragEnterEvent
        #self._ui.widget_shuttering.dropEvent = self.dropEvent
        #self._ui.widget_shuttering.dragMoveEvent = self.dragEnterEvent
        '''
        self._ui.ws = dragAbleWidget(self._ui.widget_shuttering, self)
        #self._ui.Form.__class__.dragEnterEvent = self._dragEnterEvent
        #self._ui.Form.__class__.dropEvent = self._dropEvent
        self._ui.horizontalLayout.__class__.dragEnterEvent = self.dragEnterEvent
        self._ui.horizontalLayout.__class__.dropEvent = self.dropEvent
        self._ui.widget_shuttering.setAcceptDrops(True)
        #self._ui.groupBox_shuttering.setDragEnabled(True)
        #self._ui.Form.setAcceptDrops(True)
        #self._ui.Form.setDragEnabled(True)
        '''
        self._ui.widget_shuttering.setAcceptDrops(True)

        self._ui.settings2_radioButton.hide()
        self._ui.settings3_radioButton.hide()

        self._ui.settings1_radioButton.clicked.connect(self.useSimpleSeq)
        self._ui.special_radioButton.clicked.connect(self.useSimpleSeq)
        self._ui.special2_radioButton.clicked.connect(self.useSimpleSeq)
        self._ui.settings2_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings3_radioButton.clicked.connect(self.toggleSettings)

        self._max_saved=2
        self.settings = []

        self._ui.special2_radioButton.hide()

    def dragEnterEvent(self,event):
        event.acceptProposedAction()
        event.accept()
        print "Drag Enter Event..."
        '''
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()
            '''

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            filename = str(url.encodedPath())[1:]
            print filename
            try:
                self._newSettings(filename)
            except:
                print "Unable to load settings file ", filename

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

    def toggleSettings(self):
        for i,p in enumerate(self.settings):
            radioButton = getattr(self._ui, "settings"+str(i+2)+"_radioButton")
            if radioButton.isChecked():
                self._control.updateSettings(p)
        self._ui.pushButton_startShutters.setEnabled(True)
        print "\n Shutter Sequence: "
        print p

    def useSimpleSeq(self):
        if self._ui.settings1_radioButton.isChecked() or self._ui.special_radioButton.isChecked() or self._ui.special2_radioButton.isChecked():
            self._ui.pushButton_readyShutters.setEnabled(True)
        else:
            self._ui.pushButton_readyShutters.setEnabled(False)

    def _enableShuttering(self):
        if self._ui.checkBox_enableShuttering.isChecked():
            self._ui.groupBox_shuttering401.show()
            self._ui.groupBox_shuttering488.show()
            self._ui.groupBox_shuttering561.show()
            self._ui.groupBox_shuttering642.show()
            self._ui.groupBox_shutteringLED.show()
            self.useShuttering = True
        else:
            self._ui.groupBox_shuttering401.hide()
            self._ui.groupBox_shuttering488.hide()
            self._ui.groupBox_shuttering561.hide()
            self._ui.groupBox_shuttering642.hide()
            self._ui.groupBox_shutteringLED.hide()
            self.useShuttering = False
            self._ui.pushButton_startShutters.setEnabled(False)
            self._ui.settings1_radioButton.setChecked(False)
            self._ui.settings2_radioButton.setChecked(False)
            self._ui.settings3_radioButton.setChecked(False)
            self._ui.special_radioButton.setChecked(False)
            self._ui.pushButton_readyShutters.setEnabled(False)
            
        #self._ui.pushButton_readyShutters.setEnabled(self.useShuttering)

    def _createShutterDict(self, frames):
        shutterDict = {}
        onBoxes = [self._ui.doubleSpinBox_on642,
                   self._ui.doubleSpinBox_on488,
                   self._ui.doubleSpinBox_on401,
                   self._ui.doubleSpinBox_on561,
                   self._ui.doubleSpinBox_onLED]
        offBoxes = [self._ui.doubleSpinBox_off642,
                    self._ui.doubleSpinBox_off488,
                    self._ui.doubleSpinBox_off401,
                    self._ui.doubleSpinBox_off561,
                    self._ui.doubleSpinBox_offLED]
        if self._ui.settings1_radioButton.isChecked():
            for i in range(0,len(onBoxes)):
                shutterDict[i] = {}
                on = onBoxes[i].value()
                if on>=0 and on<frames:
                    shutterDict[i]['ch_on'] = on
                else:
                    shutterDict[i]['ch_on'] = 0
                off = offBoxes[i].value()
                if off>0 and off<=frames:
                    shutterDict[i]['ch_off'] = off
                else:
                    shutterDict[i]['ch_off'] = 0
                shutterDict[i]['power'] = 1.0
            return shutterDict
        if self._ui.special_radioButton.isChecked():
            self._ui.lineEdit_shutterFrames.setText('8')
            for i in range(0,len(onBoxes)):
                shutterDict[i] = {}
                if i==0:
                    shutterDict[i]['ch_on'] = (1,5)
                    shutterDict[i]['ch_off'] = (4,8)
                elif i==2:
                    shutterDict[i]['ch_on'] = 0
                    shutterDict[i]['ch_off'] = 1
                elif i==3:
                    shutterDict[i]['ch_on'] = 4
                    shutterDict[i]['ch_off'] = 5
                else:
                    shutterDict[i]['ch_on'] = 0
                    shutterDict[i]['ch_off'] = 0
                shutterDict[i]['power'] = 1.0
            return shutterDict
        if self._ui.special2_radioButton.isChecked():
            self._ui.lineEdit_shutterFrames.setText('12')
            for i in range(0,len(onBoxes)):
                shutterDict[i] = {}
                if i==0:
                    shutterDict[i]['ch_on'] = (1,5,9)
                    shutterDict[i]['ch_off'] = (4,8,12)
                elif i==2:
                    shutterDict[i]['ch_on'] = (0,4)
                    shutterDict[i]['ch_off'] = (1,5)
                elif i==3:
                    shutterDict[i]['ch_on'] = (0,8)
                    shutterDict[i]['ch_off'] = (1,9)
                elif i==4:
                    shutterDict[i]['ch_on'] = 0.2
                    shutterDict[i]['ch_off'] = 0.8
                else:
                    shutterDict[i]['ch_on'] = 0
                    shutterDict[i]['ch_off'] = 0
                shutterDict[i]['power'] = 1.0
            return shutterDict
        return False

    def _getShutterParams(self):
        #Triggered by "Ready Shutters" button
        frames = int(self._ui.lineEdit_shutterFrames.text())
        oversampling = int(self._ui.lineEdit_oversampling.text())
        shutter_dict = self._createShutterDict(frames)
        frames = int(self._ui.lineEdit_shutterFrames.text())
        freq = int(self._ui.lineEdit_freq.text())
        print "\n Shutter sequence: "
        print shutter_dict
        useConfocal = self._ui.checkBox_confocal.isChecked()
        if useConfocal:
            #oversampling = 1000
            noCounter = True
        else:
            noCounter = False
        if shutter_dict is not False:
            wv = self._control.setShutters(5, frames, oversampling, shutter_dict, noCounter=noCounter)
            self._ui.pushButton_startShutters.setEnabled(True)
        if self._ui.checkBox_ensureDigMod.isChecked():
            self._ensureDigModCorrect(shutter_dict)
        #self.doPlot(wv)
        
    def _startShutters(self):
        useConfocal = self._ui.checkBox_confocal.isChecked()
        if useConfocal:
            noCounter = True
        else:
            noCounter = False
        if self._control.startShutters(noCounter=noCounter):
            self._ui.pushButton_startShutters.setText('STOP')
            self._ui.pushButton_startShutters.setStyleSheet("QPushButton { color: red }")
        else:
            self.stopShutters()

    def stopShutters(self):
        self._ui.pushButton_startShutters.setText('START')
        self._ui.pushButton_startShutters.setStyleSheet("QPushButton { color: black }")
        self._ui.pushButton_startShutters.setEnabled(False)
        self._control.shuttersStop()

    def loadShutterFile(self):
        shutterFile = str(QtGui.QFileDialog.getOpenFileName(self._window,
                                                            'Open File',
                                                            'D:\\'))
        if shutterFile == '' or shutterFile is None:
            return
        else:
            self._ui.labelShutterFile.setText(shutterFile)
            self._newSettings(shutterFile)                                                            

            
    def _onThroughDAQ(self, channel):
        if channel:
            state = self._ui.OnLED.isChecked()
            self._control.enableLED(state)
        else:
            state = self._ui.On561.isChecked()
            self._control.enable561(state)

    def _powerThroughDAQ(self, channel):
        if channel:
            power = self._ui.SliderLED.value()
            voltage = self._control.led_pow_func(power/100.0)
            self._control.setLED(voltage)
        else:
            power = self._ui.Slider561.value()
            voltage = self._control.sapphire_pow_func(power/100.0)
            self._control.set561(voltage)
            self._ui.Pow561.setText('Power: %.2f' % power)
            
    def _updatePower(self, index, port):
        sliderValue = self.sliders[index].value()
        power_in_mw = self._control._powerInMW(port, sliderValue/100.)
        self._control.setDiode(port, power_in_mw)
        self.powerTxts[index].setText('Power: %.2f' % sliderValue)

    def _setDigitalMod(self, index, port):
        state = self.digMod[index].isChecked()
        self._control.setDigitalMod(port, state)
        '''
        if state:
            self._control._control.lasers.setExtControl(port, "DIG")
        else:
            self._control._control.lasers.setInternalCW(port)
        '''

    def _ensureDigModCorrect(self, channelSettings):
        shouldBeDigMod = self._control.ensureDigitalModCorrect(channelSettings)
        j=2
        for i in range(len(shouldBeDigMod)):
            if shouldBeDigMod[i]:
                self.digMod[j].setChecked(True)
            j = j-1

    def _turnOnOff(self, index, port):
        state = self.onoff[index].isChecked()
        self._control.enableDiode(port, state)
        

    '''
    def doPlot(self, waveforms):
        len_wave = len(waveforms)/5
        xaxis = np.arange(0,len_wave)
        curve1 = Qwt5.QwtPlotCurve('')
        curve2 = Qwt5.QwtPlotCurve('')
        curve3 = Qwt5.QwtPlotCurve('')
        curve4 = Qwt5.QwtPlotCurve('')
        curve5 = Qwt5.QwtPlotCurve('')
        curves = [curve1, curve2, curve3, curve4, curve5]
        colors = [QtCore.Qt.red, QtCore.Qt.cyan,
                  QtCore.Qt.magenta, QtCore.Qt.yellow,
                  QtCore.Qt.black]
        pen.setStyle(QtCore.Qt.SolidLine)
        i=0
        j=0
        axisFont = QtCore.Qt.QFont('Sanserif', 6)
        self._ui.qwtPlot_shutters.setAxisFont(Qwt5.QwtPlot.xBottom, axisFont)
        self._ui.qwtPlot_shutters.setAxisFont(Qwt5.QwtPlot.yLeft, axisFont)
        for curve in curves:
            pen = QtGui.QPen(colors[j], 5)
            curve.setPen(pen)
            curve.attach(self._ui.qwtPlot_shutters)
            curve.setData(xaxis,waveforms[i:i+len_wave])
            i = i+len_wave
            j=j+1
        self._ui.qwtPlot_shutters.replot()
        for curve in curves:
            curve.detach()
    '''
