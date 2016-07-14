#!/usr/bin/python
import scipy
import inLib
from PyQt4 import QtGui, QtCore, Qwt5
from Utilities import QExtensions as qext
import numpy as np
import time

class UI(inLib.ModuleUI):
    def __init__(self, control, ui_control):
        
        inLib.ModuleUI.__init__(self, control, ui_control,
                                'modules.focusLock.focusLock_design')

        self.updaterTime = 1100
        self._ui.samples_lineEdit.setText('%i' % self._control.samples)
        self._ui.sampleRate_lineEdit.setText('%i' % self._control.sampleRate)
        self._ui.lock_lineEdit.setText('--')
        self._ui.updaterTime_lineEdit.setText('%i' % self.updaterTime)

        self._ui.samples_lineEdit.returnPressed.connect(self._setParams)
        self._ui.sampleRate_lineEdit.returnPressed.connect(self._setParams)
        self._ui.updaterTime_lineEdit.returnPressed.connect(self.setUpdaterTime)
        self._ui.prop_lineEdit.returnPressed.connect(self.setProp)
        self._ui.lock_lineEdit.returnPressed.connect(self.setLock)

        self._ui.current_pushButton.clicked.connect(self.currentLock)
        self._ui.start_pushButton.clicked.connect(self.monitorRun)
        self._ui.saveData_pushButton.clicked.connect(self.saveData)
        self._ui.saveData_pushButton.setEnabled(False)
        self._ui.clearData_pushButton.clicked.connect(self.clearData)

        self._ui.isOnFeedback_checkBox.clicked.connect(self.feedback)

        self.noPlot = False
        self.plotRealtime = True
        self._ui.radioButton_plotRealtime.setChecked(True)
        self._ui.radioButton_plotRealtime.toggled.connect(self.changePlotOption)

        #For plotting calibration data:
        grid = Qwt5.QwtPlotGrid()
        canvas = self._ui.qwtPlot.canvas().setLineWidth(2)
        pen = QtGui.QPen(QtCore.Qt.lightGray)
        grid.setPen(pen)
        grid.attach(self._ui.qwtPlot)
        self.curve = Qwt5.QwtPlotCurve('')
        self.curve.setPen(pen)
        self.curve.attach(self._ui.qwtPlot)

        self.useFeedback = False

        # A timer to update:
        self._updater = QtCore.QTimer()
        self._updater.timeout.connect(self.updateData)
        #self._updater.start(100)

        self.monitor = None
        self.errmon = None

    def _setParams(self):
        samples = int(self._ui.samples_lineEdit.text())
        sampleRate = int(self._ui.sampleRate_lineEdit.text())
        settings = {}
        settings['samples'] = samples
        settings['sampleRate'] = sampleRate
        self._control.setSettings(settings)

    def setUpdaterTime(self):
        time = int(self._ui.updaterTime_lineEdit.text())
        self.updaterTime = time

    def changePlotOption(self, state):
        if self._ui.radioButton_noPlot.isChecked():
            self.noPlot = True
        else:
            if state:
                self.plotRealtime = True
            else:
                self.plotRealtime = False

    def feedback(self):
        self.useFeedback = self._ui.isOnFeedback_checkBox.isChecked()

    def setProp(self):
        prop = float(self._ui.prop_lineEdit.text())
        self._control.setPropFeedback(prop)

    def setLock(self):
        lock = float(self._ui.lock_lineEdit.text())
        self._control.setLock(lock)

    def currentLock(self):
        current = self._control.setCurrentLock()
        if current:
            self._ui.lock_lineEdit.setText('%.3f' % current)
        else:
            self._ui.lock_lineEdit.setText('--')

    def doPlot(self, data, data_mean, means):
        if self.plotRealtime:
            xaxis = np.arange(0,len(data))
            curve1 = Qwt5.QwtPlotCurve('')
            curve2 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.black)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve1.setPen(pen)
            curve1.attach(self._ui.qwtPlot)
            curve1.setData(xaxis,data)
            pen = QtGui.QPen(QtCore.Qt.red)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve2.setPen(pen)
            curve2.attach(self._ui.qwtPlot)
            curve2.setData(np.array([xaxis[0],xaxis[-1]]),
                           np.array([data_mean, data_mean]))
            self._ui.qwtPlot.replot()
            curve1.detach()
            curve2.detach()
        else:
            xaxis = np.arange(0,len(means))
            curve1 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.black)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve1.setPen(pen)
            curve1.attach(self._ui.qwtPlot)
            curve1.setData(xaxis,means)
            self._ui.qwtPlot.replot()
            curve1.detach()

    def updateData(self):
        data,data_mean = self._control.getData(respond = self.useFeedback)
        means = self._control.getDataMeans()
        #print "Data mean: ", data_mean
        if not self.noPlot:
            self.doPlot(data, data_mean, means)
        self._ui.label_mean.setText('Mean: %.4f' % data_mean)
        self._setParams()
        self._control.startMonitor()
        
    def monitorRun(self):
        if self._ui.start_pushButton.text() == 'START':
            self._ui.start_pushButton.setText('STOP')
            self._control.startMonitor()
            time.sleep(1)
            self._updater.start(self.updaterTime)
            self._ui.saveData_pushButton.setEnabled(False)
        else:
            self._ui.start_pushButton.setText('START')
            self._updater.stop()
            self._control.resetMonitor()
            self._ui.saveData_pushButton.setEnabled(True)

    def monitorDone(self):
        self._updater.stop()
        self._ui.start_pushButton.setText('STOP')
        print 'Done...'

    def saveData(self):
        self._control.saveData()

    def clearData(self):
        self._control.reInit()

    def shutDown(self):
        self._updater.stop()
        self._control.stopMonitor()


