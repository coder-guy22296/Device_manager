#!/usr/bin/python
import scipy
import inLib
from PyQt4 import QtGui, QtCore, Qwt5
from Utilities import QExtensions as qext
import numpy as np
import time

def clickable(widget):
    class Filter(QtCore.QObject):
        clicked = QtCore.pyqtSignal(int,int)
        def eventFilter(self,obj,event):
            if obj == widget:
                if event.type() == QtCore.QEvent.MouseButtonRelease:
                    self.clicked.emit(event.x(),event.y())
                    return True
            return False
    filter=Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked

class UI(inLib.ModuleUI):
    def __init__(self, control, ui_control):
        
        inLib.ModuleUI.__init__(self, control, ui_control,
                                'modules.imageBasedFocusLock.imageBasedFocusLock_design')

        self._autoscale = False
        self._pixmap = None
        self._cmap = None
        self._vmin = 0.0
        self._vmax = 255.0
        self._previewOn = False

        self.imageXSize = 100
        self.imageYSize = 100
        self.imageX0 = 40
        self.imageY0 = 30

        self._ui.exposuretime_lineEdit.setText('%.3f' % self._control.exposure_time)
        self._ui.monitorDelay_lineEdit.setText('%.3f' % self._control.monitorDelay)
        self._ui.lineEdit_xshape.setText("%i" % self.imageXSize)
        self._ui.lineEdit_yshape.setText("%i" % self.imageYSize)
        self._ui.lineEdit_x0.setText("%i" % self.imageX0)
        self._ui.lineEdit_y0.setText("%i" % self.imageY0)
        
        self._ui.propx_lineEdit.setText('%.4f' % self._control.micronperpixel_x)
        self._ui.propy_lineEdit.setText('%.4f' % self._control.micronperpixel_y)
        self._ui.xythresh_lineEdit.setText('%.4f' % self._control.xythreshold_pixels)

        self._ui.reflectionMin_lineEdit.setText('%i' % self._control.reflectionMin)

        self._ui.propz_lineEdit.setText('%.4f' % self._control.conversion)
        self._ui.intz_lineEdit.setText('%.4f' % self._control.integrated_conv)
        self._ui.rolling_lineEdit.setText('%i' % self._control.rolling)
        self._ui.threshold_lineEdit.setText('%.4f' % self._control.threshold_diff)
        self._ui.zoff_lineEdit.setText('%.4f' % self._control.zoffset)

        self._ui.target_lineEdit.setText('%.3f' % self._control.target)
        self._ui.below_lineEdit.setText('%.3f' % self._control.below)
        self._ui.above_lineEdit.setText('%.3f' % self._control.above)

        self._ui.maybeTrouble_lineEdit.setText("%i" % self._control.maybeTroubleThreshold)

        self._ui.zsep_doubleSpinBox.setValue(0.7)
        self._ui.zsep_doubleSpinBox.setSingleStep(0.1)

        self._ui.target_lineEdit.returnPressed.connect(self.setTarget)
        self._ui.below_lineEdit.returnPressed.connect(self.setBelow)
        self._ui.above_lineEdit.returnPressed.connect(self.setAbove)

        self._ui.getTargets_pushButton.clicked.connect(self.getTargets)
        self._ui.propx_lineEdit.returnPressed.connect(self.propxyChange)
        self._ui.propy_lineEdit.returnPressed.connect(self.propxyChange)
        self._ui.xythresh_lineEdit.returnPressed.connect(self.xythreshChange)
        self._ui.propz_lineEdit.returnPressed.connect(self.propzChange)
        self._ui.intz_lineEdit.returnPressed.connect(self.intzChange)
        self._ui.threshold_lineEdit.returnPressed.connect(self.zthreshChange)
        self._ui.rolling_lineEdit.returnPressed.connect(self.rollingChange)
        self._ui.zoff_lineEdit.returnPressed.connect(self.zoffChange)
        self._ui.exposuretime_lineEdit.returnPressed.connect(self.expTimeChange)
        self._ui.calMove_lineEdit.returnPressed.connect(self.move_calibration)
        self._ui.reflectionMin_lineEdit.returnPressed.connect(self.changeReflectionMin)
        clickable(self._ui.labelDisplay).connect(self._mousePressEvent)

        self._ui.newTarget_pushButton.clicked.connect(self.newTarget)
        self._ui.targetToCurrent_pushButton.clicked.connect(self.setTargetToCurrent)

        self._ui.zsep_doubleSpinBox.valueChanged.connect(self.zsepChange)

        self._ui.noz_checkBox.stateChanged.connect(self.useZ)
        self._ui.dryrun_checkBox.stateChanged.connect(self.dryrun)
        self._ui.checkBox_noxy.stateChanged.connect(self.useXY)

        self._ui.checkBox_useMarz.setChecked(self._control.use_marz)
        self._ui.checkBox_useMarz.stateChanged.connect(self.useMarz)
        self._ui.checkBox_saveAllIms.setChecked(self._control.saveAllImages)
        self._ui.checkBox_saveAllIms.stateChanged.connect(self.saveAllIms)
        self._ui.checkBox_rollingXY.setChecked(self._control.rollingAvXY)
        self._ui.checkBox_rollingXY.stateChanged.connect(self.setRollingXY)

        self._ui.useReflection_checkBox.stateChanged.connect(self.checkedUseReflection)

        self._ui.monitorDelay_lineEdit.returnPressed.connect(self.monitorDelayChange)

        self._ui.pushButton_preview.clicked.connect(self._imagePreviewer)
        #self._ui.shutdown_pushButton.clicked.connect(self.terminate)
        self._ui.fl_pushButton.clicked.connect(self.monitorRunMP)
        self._ui.reset_pushButton.clicked.connect(self.resetFL)

        self._ui.stop_pushButton.clicked.connect(self.stopFL)

        self._ui.upperLim_lineEdit.returnPressed.connect(self._upperLim)
        self._ui.lowerLim_lineEdit.returnPressed.connect(self._lowerLim)
        self._ui.verticalSlider.valueChanged.connect(self.sliderChange)
        self._ui.setSlider_pushButton.clicked.connect(self.sliderLock)
        self._ui.verticalSlider.setMinimum(0)
        self._ui.verticalSlider.setMaximum(76800)

        self._ui.scaleValue_doubleSpinBox.setValue(0.5)

        #For plotting calibration data:
        grid = Qwt5.QwtPlotGrid()
        canvas = self._ui.qwtPlot.canvas().setLineWidth(2)
        pen = QtGui.QPen(QtCore.Qt.lightGray)
        grid.setPen(pen)
        grid.attach(self._ui.qwtPlot)
        self.curve = Qwt5.QwtPlotCurve('')
        self.curve.setPen(pen)
        self.curve.attach(self._ui.qwtPlot)

        self._ui.qwtPlot.enableAxis(1)

        self._ui.image2size_pushButton.clicked.connect(self._changeROI2)
        self._ui.use2ndImage_checkBox.stateChanged.connect(self.use2ndImage)
        self.useTwoImages = False
        self._ui.rot90_checkBox.stateChanged.connect(self._rot90)

        self._ui.maybeTrouble_lineEdit.returnPressed.connect(self._maybeTroubleThresh)

        self._ui.pushButton_flatfield.clicked.connect(self._flatfield)

        self._ui.findCOM_checkBox.stateChanged.connect(self.setFindCOM)
        self._ui.resetCOM_pushButton.clicked.connect(self.resetCOM)


        # A timer to update the image:
        self._updater = QtCore.QTimer()
        self._updater.timeout.connect(self._updateImage)
        #self._updater.start(100)

        self.monitor = None
        self.errmon = None

        self.x=0
        self.y=0

    def _changeROI(self):
        xshape = int(self._ui.lineEdit_xshape.text())
        yshape = int(self._ui.lineEdit_yshape.text())
        self._control.setXYDim(xshape,yshape)
        x0 = int(self._ui.lineEdit_x0.text())
        y0 = int(self._ui.lineEdit_y0.text())
        self._control.setX0Y0(x0,y0)

    def _changeROI2(self):
        xshape = int(self._ui.lineEdit_xshape_2.text())
        yshape = int(self._ui.lineEdit_yshape_2.text())
        x0 = int(self._ui.lineEdit_x0_2.text())
        y0 = int(self._ui.lineEdit_y0_2.text())
        self._control.setSecondImageSize(xshape,yshape,
                                         x0,y0)

    def _flatfield(self):
        if self._ui.pushButton_flatfield.text() == 'Flat field':
            self._control.flatField()
            self._ui.pushButton_flatfield.setText('End flat field')
        else:
            self._control.stopFlatField()
            self._ui.pushButton_flatfield.setText('Flat field')

    def _rot90(self):
        state = self._ui.rot90_checkBox.isChecked()
        self._control.setRotate90(state)

    def setFindCOM(self):
        state = self._ui.findCOM_checkBox.isChecked()
        self._control.setFindCOM(state)

    def resetCOM(self):
        self._control.resetCOMs()

    def _mousePressEvent(self,x,y):
        self.x = int(x)
        self.y = int(y)

    def _imagePreviewer(self):
        if not self._previewOn:
            self._changeROI()
            self.expTimeChange()
            self._control.readyCam()
            self._updater.start(100)
            self._previewOn = not self._previewOn
            self._ui.pushButton_preview.setText('Preview Stop')
        else:
            self._updater.stop()
            self._previewOn = not self._previewOn
            self._ui.pushButton_preview.setText('Preview')
                    

    def _updateImage(self):
        if self.useTwoImages:
            np_image, np_image2 = self._control.getImage(secondImage=True)
        else:
            np_image = self._control.getImage()
            np_image2 = None
        if np_image is not None:
            np_min = np_image.min()
            np_max = np_image.max()
            if self._autoscale:
                self._vmin = np_min
                self._vmax = np_max
            else:
                if np_max > 1000:
                    self._vmax = 65500
            self._ui.label_minmaxpixel.setText("Min, max pixel value: %i, %i" % (np_min, np_max))
            qt_image = qext.numpy_to_qimage8(np_image, self._vmin, self._vmax, self._cmap)
            self._pixmap = QtGui.QPixmap.fromImage(qt_image)
            self._pixmap = self._pixmap.scaled(384, 384, QtCore.Qt.KeepAspectRatio)
            self._ui.labelDisplay.setPixmap(self._pixmap)
            x = self.x*(np_image.shape[0]/384.)
            y = self.y*(np_image.shape[1]/384.)
            if (x<np_image.shape[0]) and (y<np_image.shape[1]): 
                pixelValue = np_image[y,x]
                self._ui.label_xyValue.setText("x,y: %i,%i; value: %.1f" % (x,y,pixelValue))
            else:
                self._ui.label_xyValue.setText("x,y: --,--; value: --")
            if self._ui.plotColumn_checkBox.isChecked() and (self.monitor is None) and (not self._ui.findCOM_checkBox.isChecked()):
                colSum = np_image.sum(axis=0)
                xaxis = np.arange(0,len(colSum))
                curve1 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.black)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve1.setPen(pen)
                curve1.attach(self._ui.qwtPlot)
                curve1.setData(xaxis,colSum)
                
                
                xc,xmodel = self._control.justFitPeak(colSum)
                curve2 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.blue)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve2.setPen(pen)
                curve2.attach(self._ui.qwtPlot)
                curve2.setData(xaxis,xmodel)
                self._ui.currentsig_label.setText("%.4f" % xc)
                
                self._ui.qwtPlot.replot()
                curve1.detach()
                curve2.detach()
            elif self._ui.plotColumn_checkBox.isChecked() and (self.monitor is None) and self._ui.findCOM_checkBox.isChecked():
                coms = self._control.getCOMs()
                xaxis = np.arange(0,len(coms))
                curve1 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.black)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve1.setPen(pen)
                curve1.attach(self._ui.qwtPlot)
                curve1.setData(xaxis,coms)
                self._ui.qwtPlot.replot()
                curve1.detach()
        if np_image2 is not None:
            np_min = np_image2.min()
            np_max = np_image2.max()
            if self._autoscale:
                self._vmin = np_min
                self._vmax = np_max
            qt_image = qext.numpy_to_qimage8(np_image2, self._vmin, self._vmax, self._cmap)
            self._pixmap = QtGui.QPixmap.fromImage(qt_image)
            self._pixmap = self._pixmap.scaled(256,256, QtCore.Qt.KeepAspectRatio)
            self._ui.labelDisplay2.setPixmap(self._pixmap)


    def resetFL(self):
        #print "Before reInit: len of ts: ", len(self._control.target_signal)
        self._control = self._control.reInit()
        #print "After reInit: len of ts: ", len(self._control.target_signal)
        if self.errmon is not None:
            self.errmon.stop()
            self.errmon.wait()
            self.errmon.quit()
        if self.monitor is not None:
            self.monitor.stop()
            self.monitor.wait()
        self.monitor = None
        self.errmon = None

    def newTarget(self):
        isRelative = self._ui.relative_checkBox.isChecked()
        newTarg = float(self._ui.newTarget_lineEdit.text())
        zpos = self._control.acquireNewTarget(newTarg, relative=isRelative)
        self._ui.target_lineEdit.setText('%.3f' % zpos)
        self.changeTarget()

    def _upperLim(self):
        upLim = float(self._ui.upperLim_lineEdit.text())*100.0
        self._ui.verticalSlider.setMaximum(int(upLim))

    def _lowerLim(self):
        lowLim = float(self._ui.lowerLim_lineEdit.text())*100.0
        self._ui.verticalSlider.setMinimum(int(lowLim))

    def sliderLock(self):
        '''
        lock = self._ui.verticalSlider.value()/100.0
        self._control.changeReflectionTarget(lock)
        self._ui.targetsig_label.setText('%.3f' % lock)
        '''
        lock = float(self._ui.newFocusValue_lineEdit.text())
        scale = self._ui.scaleValue_doubleSpinBox.value()
        if self.monitor is not None:
            self.monitor.setNewFL(lock,scale=scale)

    def _maybeTroubleThresh(self):
        threshold = int(self._ui.maybeTrouble_lineEdit.text())
        self._control.changeMaybeTroubleThreshold(threshold)

    def sliderChange(self):
        sliderVal = self._ui.verticalSlider.value()/100.0
        self._ui.sliderVal_label.setText('%.3f' % sliderVal)

    def move_calibration(self):
        dist = float(self._ui.calMove_lineEdit.text())
        repeats = int(self._ui.calRepeats_lineEdit.text())
        if self._ui.radioButton_x.isChecked():
            xyz = 'x'
        elif self._ui.radioButton_y.isChecked():
            xyz = 'y'
        else:
            xyz = 'z'
        useRefl = self._ui.useReflection_checkBox.isChecked()
        sigs = self._control.move_calibration(xyz, dist, repeats, reflectionBased = useRefl)
        np.save('signals_'+xyz+'cal_'+str(dist)+'.npy',sigs)
        self._ui.plotColumn_checkBox.setChecked(False)
        self.doPlot(sigs,xyz,distance=dist)
        self._ui.zoff_lineEdit.setText("%.7f" % self._control.getZOffset())

    def doPlot(self, sigs, xyz, distance=1.0):
        xaxis = np.arange(0,sigs.shape[0])
        up = np.arange(0,sigs.shape[0],2)
        down = np.arange(1,sigs.shape[0],2)
        curve1 = Qwt5.QwtPlotCurve('')
        pen = QtGui.QPen(QtCore.Qt.black)
        pen.setStyle(QtCore.Qt.SolidLine)
        curve1.setPen(pen)
        curve1.attach(self._ui.qwtPlot)
        curve2 = Qwt5.QwtPlotCurve('')
        pen = QtGui.QPen(QtCore.Qt.red)
        pen.setStyle(QtCore.Qt.DashLine)
        curve2.setPen(pen)
        curve2.attach(self._ui.qwtPlot)
        curve3 = Qwt5.QwtPlotCurve('')
        pen = QtGui.QPen(QtCore.Qt.red)
        pen.setStyle(QtCore.Qt.DashLine)
        curve3.setPen(pen)
        curve3.attach(self._ui.qwtPlot)
        if sigs.ndim == 1:
            f_sig = sigs[:]
        else:
            if xyz == 'z':
                f_sig = sigs[:,1]/sigs[:,0]
            elif xyz=='x':
                f_sig = sigs[:,2]
            elif xyz=='y':
                f_sig = sigs[:,3]
        ups = f_sig[up].mean()
        downs = f_sig[down].mean()
        diff = abs(ups-downs)
        self._ui.label_calibration.setText('Calibration: %.4f' % diff)
        self._ui.label_calResult.setText(' %.4f microns/sig' % (distance/diff))
        self._ui.propz_lineEdit.setText("%.5f" % ((distance/diff)*-0.15))
        self._ui.intz_lineEdit.setText("%.5f" % ((distance/diff)*-0.015))
        self._control.setZProp((distance/diff)*-0.15)
        self._control.setZInt((distance/diff)*-0.015)
        curve1.setData(xaxis,f_sig)
        curve2.setData(np.array([xaxis[0],xaxis[-1]]),np.array([ups, ups]))
        curve3.setData(np.array([xaxis[0],xaxis[-1]]),np.array([downs,downs]))
        self._ui.qwtPlot.replot()
        curve1.detach()
        curve2.detach()
        curve3.detach()
        

    def expTimeChange(self):
        exposure_time = float(self._ui.exposuretime_lineEdit.text())
        self._control.setExposureTime(exposure_time)

    '''
    def captDelayChange(self):
        delay = float(self._ui.captdelay_lineEdit.text())
        self._control.changeCamCaptDelay(delay)
    '''

    def use2ndImage(self):
        self.useTwoImages = self._ui.use2ndImage_checkBox.isChecked()

    def monitorDelayChange(self):
        monDelay = float(self._ui.monitorDelay_lineEdit.text())
        self._control.changeMonitorDelay(monDelay)

    def dryrun(self):
        self._control._setDryRun(self._ui.dryrun_checkBox.isChecked())

    def useXY(self):
        self._control.setNoXY(self._ui.checkBox_noxy.isChecked())

    def useMarz(self):
        self._control.setMarz(self._ui.checkBox_useMarz.isChecked())

    def saveAllIms(self):
        self._control.setSaveAll(self._ui.checkBox_saveAllIms.isChecked())

    def setRollingXY(self):
        self._control.setRollingXY(self._ui.checkBox_rollingXY.isChecked())

    def zsepChange(self,value):
        below = self._control.target - value
        above = self._control.target + value
        self._control.setTargets(above=above, below=below)
        self._ui.below_lineEdit.setText('%.3f' % below)
        self._ui.above_lineEdit.setText('%.3f' % above)

    def changeReflectionMin(self):
        refMin = int(self._ui.reflectionMin_lineEdit.text())
        self._control.setReflectionMin(refMin)

    def changeTarget(self):
        targetPos = float(self._ui.target_lineEdit.text())
        self._control.setTargets(target=targetPos)
        
    def useZ(self,val):
        self._control.no_z = self._ui.noz_checkBox.isChecked()
        
    def getTargets(self):
        self.setTarget()
        self.setAbove()
        self.setBelow()
        self.gtThread = getTargetsThread(self._control)
        self._window.connect(self.gtThread, QtCore.SIGNAL('gotTargets'), self.gotTargets)
        self.gtThread.start()
        self._ui.getTargets_pushButton.setEnabled(False)
        #self._control.acquireTargets()
        #self._control.fftReferences()

    def gotTargets(self):
        self._ui.calMove_lineEdit.setEnabled(True)
        self._ui.fl_pushButton.setEnabled(True)
        self._ui.getTargets_pushButton.setEnabled(True)
        self.gtThread.wait(100)

    def checkedUseReflection(self):
        if self._ui.useReflection_checkBox.isChecked():
            self._ui.fl_pushButton.setEnabled(True)
            self._ui.calMove_lineEdit.setEnabled(True)
            self._ui.tab_xy.setEnabled(False)
            [propz, intz, thresh, rolling] = self._control.setReflectionDefaultParams()
            self._ui.propz_lineEdit.setText("%.5f" % propz)
            self._ui.intz_lineEdit.setText("%.6f" % intz)
            self._ui.threshold_lineEdit.setText("%.3f" % thresh)
            self._ui.rolling_lineEdit.setText("%i" % rolling)
        else:
            self._ui.fl_pushButton.setEnabled(False)
            self._ui.calMove_lineEdit.setEnabled(False)
            self._ui.tab_xy.setEnabled(True)

    def setTarget(self):
        target = float(self._ui.target_lineEdit.text())
        self._control.setTargets(target=target)

    def setAbove(self):
        above = float(self._ui.above_lineEdit.text())
        self._control.setTargets(above=above)

    def setBelow(self):
        below = float(self._ui.below_lineEdit.text())
        self._control.setTargets(below=below)

    def setTargetToCurrent(self):
        zsep = float(self._ui.zsep_doubleSpinBox.value())
        zpos = self._control.setTargetToCurrent(zsep)
        self._ui.target_lineEdit.setText('%.4f' % zpos)
        above = zpos+zsep
        below = zpos-zsep
        self._ui.above_lineEdit.setText('%.3f' % above)
        self._ui.below_lineEdit.setText('%.3f' % below)

    def propxyChange(self):
        xprop = float(self._ui.propx_lineEdit.text())
        yprop = float(self._ui.propy_lineEdit.text())
        self._control.setXYProp(xprop,yprop)

    def xythreshChange(self):
        xythresh = float(self._ui.xythresh_lineEdit.text())
        self._control.setXYThresh(xythresh)

    def propzChange(self):
        zprop = float(self._ui.propz_lineEdit.text())
        self._control.setZProp(zprop)

    def intzChange(self):
        intz = float(self._ui.intz_lineEdit.text())
        self._control.setZInt(intz)

    def zthreshChange(self):
        zthresh = float(self._ui.threshold_lineEdit.text())
        self._control.setZThresh(zthresh)

    def rollingChange(self):
        rolling = int(self._ui.rolling_lineEdit.text())
        self._control.setRolling(rolling)

    def zoffChange(self):
        zoff = float(self._ui.zoff_lineEdit.text())
        self._control.setZOffset(zoff)

    def monitorRun(self):
        if self._ui.fl_pushButton.text() == 'Focus Lock Start':
            self._control.setQuit(False)
            self._ui.fl_pushButton.setText('Focus Lock Stop')
            num_frames = int(self._ui.numFrames_lineEdit.text())
            self._control.setFrames(num_frames)
            if num_frames>0:
                self.monitor = Monitor(self._control)
                self.errmon = ErrorMonitor(self._control, self.monitor)
                self._window.connect(self.errmon, QtCore.SIGNAL('zStageError'), self.updateZStageError)
                self._window.connect(self.monitor, QtCore.SIGNAL('monitorDone'), self.monitorDone)
                self._window.connect(self.errmon, QtCore.SIGNAL('plotErrorSignal'), self.plotErrorSignal)
                self.monitor.start()
                self.errmon.start()
        else:
            self._control.setQuit(True)
            time.sleep(0.2)
            self.monitor.stop()
            self.monitor.wait(1000)
            self.errmon.stop()
            self.errmon.wait(1000)
            self.monitor = None
            self.errmon = None
            self._ui.fl_pushButton.setText('Focus Lock Start')

    def monitorRunMP(self):
        if self._ui.fl_pushButton.text() == 'Focus Lock Start':
            if self._ui.toMove_checkBox.isChecked():
                timesToMoveStr = self._ui.timesToMove_lineEdit.text()
                positionsToMoveStr = self._ui.posToMove_lineEdit.text()
                timesToMove = []
                positionsToMove = []
                tempTimes = timesToMoveStr.split(',')
                tempPos = positionsToMoveStr.split(',')
                for k in range(0,len(tempTimes)):
                    timesToMove.append(int(tempTimes[k]))
                    positionsToMove.append(float(tempPos[k]))
                tAndPs = [timesToMove,positionsToMove]
            else:
                tAndPs = None
            self._control.setQuit(False)
            self._ui.fl_pushButton.setText('Focus Lock Stop')
            num_frames = int(self._ui.numFrames_lineEdit.text())
            self._control.setFrames(num_frames)
            if num_frames>0:
                if self._ui.useReflection_checkBox.isChecked():
                    self.monitor = MonitorReflection(self._control, num_frames, timesAndPositions=tAndPs)
                    self._window.connect(self.monitor, QtCore.SIGNAL('zStageError(float,float)'),
                                         self.updateZReflectionError)
                    self._window.connect(self.monitor, QtCore.SIGNAL('monitorDone'), self.monitorMPDone)
                    self._window.connect(self.monitor, QtCore.SIGNAL('plotModels'), self.plotModels)
                else:
                    self.monitor = MonitorMP(self._control, num_frames)
                    self._window.connect(self.monitor, QtCore.SIGNAL('zStageError(float,float,float,float,float)'),
                                         self.updateZStageError)
                    self._window.connect(self.monitor, QtCore.SIGNAL('monitorDone'), self.monitorMPDone)
                    self._window.connect(self.monitor, QtCore.SIGNAL('plotErrorSignal'), self.plotErrorSignal)
                self.monitor.start()
        else:
            self._control.setQuit(True)
            self.monitor.setInactive()
            time.sleep(0.2)
            self.monitor.stop()
            self.monitor.wait(1000)
            self.monitor = None
            self._ui.fl_pushButton.setText('Focus Lock Start')

    def reflectionBasedMonitor(self):
        if self._ui.fl_pushButton.text() == 'Focus Lock Start':
            self._control.setQuit(False)
            self._ui.fl_pushButton.setText('Focus Lock Stop')
            num_frames = int(self._ui.numFrames_lineEdit.text())
            self._control.setFrames(num_frames)
            if num_frames>0:
                self.monitor = MonitorReflection(self._control, num_frames)
                self._window.connect(self.monitor, QtCore.SIGNAL('zStageError(float,float)'),
                                     self.updateZReflectionError)
                self._window.connect(self.monitor, QtCore.SIGNAL('monitorDone'), self.monitorMPDone)
                self._window.connect(self.monitor, QtCore.SIGNAL('plotModels'), self.plotModels)
                self.monitor.start()
        else:
            self._control.setQuit(True)
            self.monitor.setInactive()
            time.sleep(0.2)
            self.monitor.stop()
            self.monitor.wait(1000)
            self.monitor = None
            self._ui.fl_pushButton.setText('Focus Lock Start')
            self.resetFL()

    def monitorDone(self):
        self._ui.fl_pushButton.setText('Focus Lock Start')
        print 'Done...'

    def monitorMPDone(self):
        self._ui.fl_pushButton.setText('Focus Lock Start')
        print 'Done...'

    '''
    def updateZStageError(self):
        #print 'xdrift len: ', len(self._control.xdrift)
        if len(self._control.stage_zs)>0 and len(self._control.xdrift)>0:
            zstage = self._control.stage_zs[-1]
            self._ui.zstage_label.setText('%.4f' % zstage)
            xd,yd = self._control.xdrift[-1], self._control.ydrift[-1]
            self._ui.xyd_label.setText('%.4f, %.4f' % (xd,yd))
        beginFrames = self._control.initialFrames + self._control.unactive_until + 1
        if len(self._control.stage_zs)>beginFrames:
            ts0 = self._control.target_signal[self._control.unactive_until:beginFrames]
            ds0 = self._control.defocus_signal[self._control.unactive_until:beginFrames]
            targetsig = np.mean(ds0)/np.mean(ts0)
            self._ui.targetsig_label.setText('%.4e' % targetsig)
            ts = self._control.target_signal[-1]
            ds = self._control.defocus_signal[-1]
            sig = float(ds)/float(ts)
            self._ui.currentsig_label.setText('%.4e' % sig)
            zm = self._control.movedz[-1]
            self._ui.ztravel_label.setText('%.4f' % zm)
    '''

    def updateZStageError(self, zstage, xd, yd, sig0, sig):
        self._ui.zstage_label.setText('%.4f' % zstage)
        self._ui.xyd_label.setText('%.4f, %.4f' % (xd,yd))
        if sig0==-1:
            self._ui.targetsig_label.setText('--')
        else:
            self._ui.targetsig_label.setText('%.4e' % sig0)
        self._ui.currentsig_label.setText('%.4e' % sig)

    def updateZReflectionError(self, xc, zstage):
        beginFrames = self._control.initialFrames + self._control.unactive_until + 1
        if len(self._control.stage_zs)<beginFrames:
            self._ui.targetsig_label.setText('--')
        elif (len(self._control.stage_zs)>beginFrames) and (str(self._ui.targetsig_label.text())=='--'):
            ts0 = np.mean(self._control.target_signal[self._control.unactive_until:beginFrames])
            self._ui.targetsig_label.setText('%.3f' % ts0)
        self._ui.currentsig_label.setText('%.3f' % xc)
        self._ui.zstage_label.setText('%.4f' % zstage)


    def plotErrorSignal(self):
        beginFrames = self._control.initialFrames + self._control.unactive_until + 1
        if len(self._control.stage_zs)>beginFrames:
            ts0 = self._control.target_signal[self._control.unactive_until:beginFrames]
            ds0 = self._control.defocus_signal[self._control.unactive_until:beginFrames]
            targetsig = np.mean(ds0)/np.mean(ts0)
            ts = self._control.target_signal
            ds = self._control.defocus_signal
            sig = np.array(ds)/np.array(ts)
            xaxis = np.arange(0,len(sig))
            curve1 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.black)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve1.setPen(pen)
            curve1.attach(self._ui.qwtPlot)
            curve2 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.red)
            pen.setStyle(QtCore.Qt.DashLine)
            curve2.setPen(pen)
            curve2.attach(self._ui.qwtPlot)
            curve3 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.blue)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve3.setPen(pen)
            curve3.setYAxis(1)
            curve3.attach(self._ui.qwtPlot)
            curve3.setData(xaxis,self._control.stage_zs)
            curve1.setData(xaxis,sig)
            curve2.setData(np.array([xaxis[0],xaxis[-1]]), np.array([targetsig,targetsig]))
            self._ui.qwtPlot.replot()
            curve1.detach()
            curve2.detach()
            curve3.detach()

    def plotModels(self):
        if self._ui.plotFit_checkBox.isChecked():
            xmodel,sumx = self._control.returnModels()
            if len(xmodel)>0:
                xaxis = np.arange(0,len(xmodel))
                curve1 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.black)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve1.setPen(pen)
                curve1.attach(self._ui.qwtPlot)
                curve1.setData(xaxis,xmodel)
                curve2 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.blue)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve2.setPen(pen)
                curve2.attach(self._ui.qwtPlot)
                curve2.setData(xaxis,sumx)
                self._ui.qwtPlot.replot()
                curve1.detach()
                curve2.detach()
        else:
            beginFrames = self._control.initialFrames + self._control.unactive_until + 1
            if len(self._control.stage_zs)>beginFrames:
                times1 = np.array(self._control.times) - self._control.times[0]
                times2 = np.array(self._control.times2) - self._control.times2[0]
                ts0 = self._control.target_signal[self._control.unactive_until:beginFrames]
                targetsig = np.mean(ts0)
                sig = np.array(self._control.target_signal)
                #xaxis = np.arange(0,len(sig))
                curve1 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.black)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve1.setPen(pen)
                curve1.attach(self._ui.qwtPlot)
                curve2 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.red)
                pen.setStyle(QtCore.Qt.DashLine)
                curve2.setPen(pen)
                curve2.attach(self._ui.qwtPlot)
                curve3 = Qwt5.QwtPlotCurve('')
                pen = QtGui.QPen(QtCore.Qt.blue)
                pen.setStyle(QtCore.Qt.SolidLine)
                curve3.setPen(pen)
                curve3.setYAxis(1)
                curve3.attach(self._ui.qwtPlot)
                curve3.setData(times1,self._control.stage_zs)
                curve1.setData(times2,sig)
                curve2.setData(np.array([times2[0],times2[-1]]), np.array([targetsig,targetsig]))
                self._ui.qwtPlot.replot()
                curve1.detach()
                curve2.detach()
                curve3.detach()
            

    def shutDown(self):
        self._updater.stop()
        if self.monitor is not None:
            self.monitor.stop()
            self.monitor.wait(2000)

    def stopFL(self):
        try:
            self._control.setQuit(True)
        except:
            print "Busy..."


class Monitor(QtCore.QThread):
    def __init__(self, focuslock):
        QtCore.QThread.__init__(self)
        self.focuslock = focuslock
        self.focuslock.active = True
        self.active=True
        print "Starting monitor..."

    def run(self):
        self.focuslock.monitor_limitedFrames_fast()
        self.emit(QtCore.SIGNAL('monitorDone'))
        self.active=False

    def stop(self):
        print "Pressed stop..."
        self.focuslock.setQuit(True)
        self.focuslock.active=False
        self.active=False

class getTargetsThread(QtCore.QThread):
    def __init__(self, focuslock):
        QtCore.QThread.__init__(self)
        self.focuslock = focuslock

    def run(self):
        self.focuslock.acquireTargets()
        self.focuslock.fftReferences()
        self.emit(QtCore.SIGNAL('gotTargets'))

    def stop(self):
        pass
        

class MonitorMP(QtCore.QThread):
    def __init__(self, focuslock, num_frames, timesAndPositions=None):
        QtCore.QThread.__init__(self)
        self.focuslock = focuslock
        self.focuslock.active = True
        self.active=True
        self.numFrames = num_frames
        print "Starting monitor..."

    def run(self):
        j=0
        while (j<self.numFrames) and self.active:
            stagezs, xd, yd, sig0, sigNow = self.focuslock.monitor_limitedFrames_MP(j)
            j = j+1
            self.active = not self.focuslock.getQuitState()
            self.emit(QtCore.SIGNAL('zStageError(float,float,float,float,float)'), stagezs, xd, yd, sig0, sigNow)
            if np.mod(j,10):
                self.emit(QtCore.SIGNAL('plotErrorSignal'))
        self.emit(QtCore.SIGNAL('monitorDone'))
        self.focuslock.finishUpFocusLock()
        self.active=False

    def setInactive(self):
        self.active=False

    def stop(self):
        self.focuslock.setQuit(True)
        self.focuslock.active=False
        self.active=False

class MonitorReflection(QtCore.QThread):
    def __init__(self, focuslock, num_frames, timesAndPositions=None):
        QtCore.QThread.__init__(self)
        self.timesAndPositions = timesAndPositions
        self.focuslock = focuslock
        self.focuslock.active = True
        self.active=True
        self.numFrames = num_frames
        print "Starting monitor..."

    def run(self):
        j=0
        timeAndPosIndex = 0
        while (j<self.numFrames) and self.active:
            xc,modelx,zstage,timenow = self.focuslock.monitor_limitedFrames_Reflection(j, timesAndPositions=self.timesAndPositions)
            j = j+1
            self.active = not self.focuslock.getQuitState()
            self.emit(QtCore.SIGNAL('zStageError(float,float)'), xc, zstage)
            if np.mod(j,10):
                self.emit(QtCore.SIGNAL('plotModels'))
            if self.timesAndPositions is not None:
                if timenow>self.timesAndPositions[0][timeAndPosIndex]:
                    self.setNewFL(self.timesAndPositions[1][timeAndPosIndex])
                    timeAndPosIndex += 1
                    if timeAndPosIndex == len(self.timesAndPositions[0]):
                        self.timesAndPositions = None
        self.emit(QtCore.SIGNAL('monitorDone'))
        self.focuslock.finishUpFocusLock(reflection=True)
        self.active=False

    def setInactive(self):
        self.active=False

    def setNewFL(self, value, scale=1.5):
        self.focuslock.changeReflectionTarget(value, scale)
        
    def stop(self):
        self.focuslock.setQuit(True)
        self.focuslock.active=False
        self.active=False

'''
class ErrorMonitor(QtCore.QThread):
    def __init__(self, focuslock, mon):
        QtCore.QThread.__init__(self)
        self.focuslock = focuslock
        self.monitor = mon
        self.active=True
        self.iter = 0
        print "Starting error monitor..."

    def run(self):
        while (self.monitor.active and self.active):
            #print "ts len: ", len(self.focuslock.target_signal)
            #if len(self.focuslock.target_signal)>0:
            #    self.emit(QtCore.SIGNAL('zStageError'))
            time.sleep(2.)
            self.iter+=1
            if np.mod(self.iter,4)==0:
                self.emit(QtCore.SIGNAL('plotErrorSignal'))

    def stop(self):
        self.focuslock.active=False
        self.active=False
'''
        
