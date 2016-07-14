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

        self.updaterTime = 110
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

        self._ui.twochannels_checkBox.clicked.connect(self.twoChannels)

        self._ui.piezoOffset_lineEdit.setText("0.0")
        self._ui.piezoOffset_lineEdit.returnPressed.connect(self._setPiezoOffset)
        self._ui.findPiezoOffset_pushButton.clicked.connect(self._findOffset)

        self._ui.calibrate_pushButton.clicked.connect(self.calibrate)

        self._ui.upperLim_lineEdit.returnPressed.connect(self._upperLim)
        self._ui.lowerLim_lineEdit.returnPressed.connect(self._lowerLim)
        self._ui.verticalSlider.valueChanged.connect(self.sliderLock)

        self.noPlot = False
        self.plotRealtime = True
        self._ui.radioButton_plotRealtime.setChecked(True)
        self._ui.radioButton_plotRealtime.toggled.connect(self.changePlotOption)
        self._ui.radioButton_noPlot.toggled.connect(self.changePlotOption)

        self._ui.setThreshold_pushButton.clicked.connect(self._setThreshold)
        

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

        self.oneAIChannel = True

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

    def _setThreshold(self):
        threshold = float(self._ui.intThreshold_lineEdit.text())
        self._control.setIntensityThreshold(threshold)

    def _upperLim(self):
        upLim = float(self._ui.upperLim_lineEdit.text())
        self._ui.verticalSlider.setMaximum(int(upLim*1000))

    def _lowerLim(self):
        lowLim = float(self._ui.lowerLim_lineEdit.text())
        self._ui.verticalSlider.setMinimum(int(lowLim*1000))
        print "Low of %.4f" % lowLim

    def sliderLock(self):
        lock = self._ui.verticalSlider.value()/1000.0
        self._control.setLock(lock)
        self._ui.lock_lineEdit.setText("%.4f" % lock)

    def calibrate(self):
        dists1, vals1, slope1 = self._control.calibrate(0.5)
        self._ui.dist1_label.setText("Distance: %.4f" % 0.5)
        self._ui.slope1_label.setText("Slope: %.4f" % slope1)
        dists2, vals2, slope2 = self._control.calibrate(1.0)
        self._ui.dist2_label.setText("Distance: %.4f" % 1.0)
        self._ui.slope2_label.setText("Slope: %.4f" % slope2)
        dists3, vals3, slope3 = self._control.calibrate(1.5)
        self._ui.dist3_label.setText("Distance: %.4f" % 1.5)
        self._ui.slope3_label.setText("Slope: %.4f" % slope3)
        self.doCalPlot(np.vstack((dists1, dists2, dists3)),
                       np.vstack((vals1, vals2, vals3)))

    def doCalPlot(self, xdata, ydata):
        curve1 = Qwt5.QwtPlotCurve('')
        curve2 = Qwt5.QwtPlotCurve('')
        curve3 = Qwt5.QwtPlotCurve('')
        pen = QtGui.QPen(QtCore.Qt.black)
        pen.setStyle(QtCore.Qt.SolidLine)
        curve1.setPen(pen)
        curve1.attach(self._ui.qwtPlot)
        curve1.setData(xdata[0],ydata[0])
        pen = QtGui.QPen(QtCore.Qt.red)
        pen.setStyle(QtCore.Qt.SolidLine)
        curve2.setPen(pen)
        curve2.attach(self._ui.qwtPlot)
        curve2.setData(xdata[1],ydata[1])
        pen = QtGui.QPen(QtCore.Qt.blue)
        pen.setStyle(QtCore.Qt.SolidLine)
        curve3.setPen(pen)
        curve3.attach(self._ui.qwtPlot)
        curve3.setData(xdata[2],ydata[2])
        self._ui.qwtPlot.replot()
        curve1.detach()
        curve2.detach()
        curve3.detach()

    def _findOffset(self):
        offset = self._control.findOffset()
        self._ui.piezoOffset_lineEdit.setText("%.4f" % offset)

    def _setPiezoOffset(self):
        offset = float(self._ui.piezoOffset_lineEdit.text())
        self._control.setOffset(offset)

    def changePlotOption(self, state):
        if self._ui.radioButton_noPlot.isChecked():
            self.noPlot = True
        elif self._ui.radioButton_plotRealtime.isChecked():
            self.plotRealtime = True
            self.noPlot = False
        else:
            self.plotRealtime = False
            self.noPlot = False

    def feedback(self):
        self.useFeedback = self._ui.isOnFeedback_checkBox.isChecked()

    def twoChannels(self):
        state = self._ui.twochannels_checkBox.isChecked()
        self._control.useTwoChannels(state)
        self.oneAIChannel = not state

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
            self._ui.verticalSlider.setValue(int(current * 1000))
        else:
            self._ui.lock_lineEdit.setText('--')

    def doPlot(self, data, data_mean, means, means2):
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
            curve2 = Qwt5.QwtPlotCurve('')
            curve3 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.black)
            pen.setStyle(QtCore.Qt.SolidLine)
            pen2 = QtGui.QPen(QtCore.Qt.blue)
            pen2.setStyle(QtCore.Qt.SolidLine)
            pen3 = QtGui.QPen(QtCore.Qt.red)
            pen3.setStyle(QtCore.Qt.SolidLine)
            curve1.setPen(pen)
            curve1.attach(self._ui.qwtPlot)
            if not self._ui.plotOnlyPos_checkBox.isChecked():
                curve1.setData(xaxis,means)
            curve2.setPen(pen2)
            curve2.attach(self._ui.qwtPlot)
            curve3.setPen(pen3)
            curve3.attach(self._ui.qwtPlot)
            if means2[0] < 100:
                xaxis2 = np.arange(0,len(means2))
                if not self._ui.plotOnlyPos_checkBox.isChecked():
                    curve2.setData(xaxis2,means2)
                curve3.setData(xaxis2,np.array(means)/(2.0*np.array(means2)))
            self._ui.qwtPlot.replot()
            curve1.detach()
            curve2.detach()
            curve3.detach()

    def updateData(self):
        if self.oneAIChannel:
            data,data_mean = self._control.getData(respond = self.useFeedback)
        else:
            data,data_mean,data_mean2 = self._control.getData(respond = self.useFeedback)
        if self.oneAIChannel:
            means = self._control.getDataMeans()
            means2 = [100]
        else:
            means, means2 = self._control.getDataMeans()
        #print "Data mean: ", data_mean
        if not self.noPlot:
            self.doPlot(data, data_mean, means, means2)
        self._ui.label_mean.setText('Mean: %.6f' % data_mean)
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


