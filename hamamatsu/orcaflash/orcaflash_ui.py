#!/usr/bin/python

import inLib
from PyQt4 import QtCore, QtGui, Qwt5
from Utilities import QExtensions as qext
from functools import partial
import time
import numpy as np
import os

def clickable(widget):
    class Filter(QtCore.QObject):
        clicked = QtCore.pyqtSignal(int,int)
        def eventFilter(self,obj,event):
            if obj == widget:
                if event.type() == QtCore.QEvent.MouseButtonRelease:
                    buttonState = event.button()
                    self.clicked.emit(event.x(),event.y())
                    return True
            return False
    filter=Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked

class UI(inLib.DeviceUI):

    def __init__(self, control):
        inLib.DeviceUI.__init__(self, control, 'hamamatsu.orcaflash.orcaflash_design')

        self._control.beginPreview()

        self._autoscale = True
        self._autoscale90 = False
        self._pixmap = None
        self._cmap = None
        self.vmin = 0.0
        self.vmax = 65535

        self.filename = 'movie0001'

        self._previewMode = True
        
        #Connecting signals
        self._ui.checkBoxAutoscale.stateChanged.connect(self.setAutoscale)          #autoscale
        self._ui.checkBox_autoscale90.stateChanged.connect(self.setAutoscale90)     #autoscale to 90% of max
        self._ui.pushButton_stopCapture.clicked.connect(self.stopCapture)           #
        self._ui.setROI_pushButton.clicked.connect(self._setROI)                    #sets region of interest
        self._ui.resetPreview_pushButton.clicked.connect(self._resetPreview)
        self._ui.pushButton_SaveFrame.clicked.connect(self._saveImageNumpy)         #saves current image as Numpy array
        self._ui.pushButton_SetTimings.clicked.connect(self._setTimings)            #set exposure/delay times
        self._ui.pushButton_Record.clicked.connect(self._record)
        self._ui.comboBox_recordMode.currentIndexChanged.connect(self._setRecordMode)
        self._ui.lineEdit_fileName.textChanged.connect(self._setFilename)
        self._ui.spinBox_fileNameIndex.valueChanged.connect(self._setFilename)
        self._ui.pushButton_printCamProps.clicked.connect(self._printCamProps)
        self._ui.pushButton_propSet.clicked.connect(self._setProp)
        self._ui.pushButton_propGet.clicked.connect(self._getProp)
        self._ui.radioButton_2x2Binning.clicked.connect(self._setBinning)
        self._ui.radioButton_noBinning.clicked.connect(self._setBinning)
        self._ui.labelDisplay.paintEvent = self._labelDisplay_paintEvent
        clickable(self._ui.labelDisplay).connect(self._mousePressEvent)
        self._ui.pushButton_beginSnaps.clicked.connect(self._beginSnaps)
        self._ui.pushButton_takeSnap.clicked.connect(self._takeSnap)
        self._ui.pushButton_doneSnaps.clicked.connect(self._doneSnaps)
        self._ui.checkBox_lowLight.stateChanged.connect(self._lowLight)
        self._ui.checkBox_defectMode.stateChanged.connect(self._defectMode)
        self._ui.lineEdit_toDiskTime.returnPressed.connect(self._newUpdateTime)
        self._ui.triggerMode_lineEdit.returnPressed.connect(self._triggerModeUpdate)
        self._ui.pushButton_triggermodes.clicked.connect(self._triggerModeUpdate)
        self._ui.pushButton_roi128128.clicked.connect(partial(self._setROItoCommon, "128x128"))
        self._ui.pushButton_roi256256.clicked.connect(partial(self._setROItoCommon, "256x256"))
        self._ui.pushButton_roi512256.clicked.connect(partial(self._setROItoCommon, "512x256"))
        self._ui.pushButton_roi512512.clicked.connect(partial(self._setROItoCommon, "512x512"))
        
        

        #ComboBoxes for exposure (and delay though perhaps NA)
        self._ui.comboBox_expunits.addItem("s")
        self._ui.comboBox_expunits.addItem("ms")
        self._ui.comboBox_expunits.addItem("us")
        self._ui.comboBox_expunits.addItem("ns")
        self._ui.comboBox_expunits.setCurrentIndex(1)
        #self._ui.comboBox_delayunits.addItem("s")
        #self._ui.comboBox_delayunits.addItem("ms")
        #self._ui.comboBox_delayunits.addItem("us")
        #self._ui.comboBox_delayunits.addItem("ns")
        #self._ui.comboBox_delayunits.setCurrentIndex(1)

        #Setting text:
        xdim,ydim = self._control._getResolution()
        self._ui.labelResolution.setText("%i, %i" % (xdim,ydim))
        self._control._props['dimensions'] = xdim,ydim
        self._xres, self._yres = xdim,ydim
        self._ui.lineEdit_fileName.setText("movie")

        self.plotWidth = 0
        self.plotLength = 0
        self.xForPlot = 1
        self.yForPlot = 1
        self._ui.lineEdit_plotWidth.setText("%i" % self.plotWidth)
        self._ui.lineEdit_plotLength.setText("%i" % self.plotLength)
        self._ui.lineEdit_xForPlot.setText("%i" % self.xForPlot)
        self._ui.lineEdit_yForPlot.setText("%i" % self.yForPlot)
        self._ui.lineEdit_plotWidth.returnPressed.connect(partial(self._setPlot,'w'))
        self._ui.lineEdit_plotLength.returnPressed.connect(partial(self._setPlot,'l'))
        self._ui.lineEdit_xForPlot.returnPressed.connect(partial(self._setPlot,'x'))
        self._ui.lineEdit_yForPlot.returnPressed.connect(partial(self._setPlot,'y'))
        self._ui.radioButton_plotHor.setChecked(True)
        

        self._getDisplayTimings()
        self._writeROIValues()

        self.x=0
        self.y=0

        ##Settings
        self._ui.settings2_radioButton.hide()
        self._ui.settings3_radioButton.hide()
        self._ui.settings4_radioButton.hide()
        self._ui.settings5_radioButton.hide()
        self._ui.settings6_radioButton.hide()

        self.settings = [self._control.loadedSettings]
        print self.settings[0]
        settings_name = self.settings[0]['settings_filename']
        self._ui.settings1_radioButton.setText(settings_name)
        self._ui.settings1_radioButton.setChecked(True)

        self._ui.settings1_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings2_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings3_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings4_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings5_radioButton.clicked.connect(self.toggleSettings)
        self._ui.settings6_radioButton.clicked.connect(self.toggleSettings)

        self._max_saved = 6

        self._ui.pushButton_5hz.clicked.connect(partial(self.commonTiming, 5))
        self._ui.pushButton_10hz.clicked.connect(partial(self.commonTiming, 10))
        self._ui.pushButton_20hz.clicked.connect(partial(self.commonTiming, 20))
        self._ui.pushButton_60hz.clicked.connect(partial(self.commonTiming, 60))
        self._ui.pushButton_100hz.clicked.connect(partial(self.commonTiming, 100))
        self._ui.pushButton_120hz.clicked.connect(partial(self.commonTiming, 120))

        ##Handling file drops
        
        self._ui.tabWidget.__class__.dragEnterEvent = self._dragEnterEvent
        self._ui.tabWidget.__class__.dropEvent = self._dropEvent
        self._ui.tabWidget.setAcceptDrops(True)
        #self._ui.tabWidget.setDragEnabled(True)
        
        self._ui.labelDisplay.setGeometry(10,10,512,512)

        #For testing purposes. Hide if not in testing mode
        #self._ui.groupBox_testing.hide()

        # A timer to update the image:
        self._updater = QtCore.QTimer()
        self._updater.timeout.connect(self._update)
        self._updater.start(50)

        self._updateToDisk = QtCore.QTimer()
        self._updateToDisk.timeout.connect(self._fillDAX)
        self._toDiskUpdateTime = 200

        self._recordingOn = False

        self._ui.spinBox_fileNameIndex.setValue(1)

        self._ui.horizontalSlider_min.setRange(0,65535)
        self._ui.horizontalSlider_max.setRange(0,65535)
        self._ui.horizontalSlider_min.valueChanged.connect(self._setDispMin)
        self._ui.horizontalSlider_max.valueChanged.connect(self._setDispMax)

        grid = Qwt5.QwtPlotGrid()
        canvas = self._ui.qwtPlot.canvas().setLineWidth(2)
        pen = QtGui.QPen(QtCore.Qt.lightGray)
        grid.setPen(pen)
        grid.attach(self._ui.qwtPlot)
        self.curve = Qwt5.QwtPlotCurve('')
        self.curve.setPen(pen)
        self.curve.attach(self._ui.qwtPlot)
        self._ui.qwtPlot.enableAxis(1)


    def _printCamProps(self):
        self._control.printProperties()

    def _newUpdateTime(self):
        uptime = int(self._ui.lineEdit_toDiskTime.text())
        if uptime>10 and uptime<10000:
            self._toDiskUpdateTime = uptime

    def _triggerModeUpdate(self):
        #value = int(self._ui.triggerMode_lineEdit.text())
        index = self._ui.comboBox_triggermodes.currentIndex()
        value = 0x0001
        if index==0:
            value = 0x0001
        elif index==1:
            value = 0x0002
        elif index==2:
            value = 0x0004
        elif index==3:
            value = 0x0008
        elif index==4:
            value = 0x0010
        elif index==5:
            value = 0x0020
        elif index==6:
            value = 0x0040
        elif index==7:
            value = 0x0080
        elif index==8:
            value = 0x0100
        elif index==9:
            value = 0x0200
        elif index==10:
            value = 0x0400
        self._control.setTriggerMode(value)

    def _lowLight(self, state):
        if state:
            self._ui.horizontalSlider_min.setRange(0,10000)
            self._ui.horizontalSlider_max.setRange(0,10000)
        else:
            self._ui.horizontalSlider_min.setRange(0,65535)
            self._ui.horizontalSlider_max.setRange(0,65535)

    def _defectMode(self, state):
        if state:
            self._control.setDefectCorrect(2)
        else:
            self._control.setDefectCorrect(1)

    def _mousePressEvent(self,x,y):
        self.x = int(x*(self._yres/512.))
        self.y = int(y*(self._xres/512.))

    def _setDispMin(self):
        self.vmin = self._ui.horizontalSlider_min.value()
        self._ui.label_minSlider.setText("Min: %.1f" % self.vmin)

    def _setDispMax(self):
        self.vmax = self._ui.horizontalSlider_max.value()
        self._ui.label_maxSlider.setText("Max: %.1f" % self.vmax)

    def commonTiming(self, hertz):
        exp_time = 1./hertz
        self._control.setExposureTime(exp_time)
        self._getDisplayTimings()
        
    def _setProp(self):
        propID = int(self._ui.lineEdit_propID.text())
        propVal = float(self._ui.lineEdit_propVal.text())
        self._control._setProp(propID,propVal)

    def _getProp(self):
        propID = int(self._ui.lineEdit_propID.text())
        propVal = self._control._getProp(propID)
        self._ui.lineEdit_propVal.setText("%.5f" % propVal)

    def _dragEnterEvent(self,event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def _dropEvent(self, event):
        for url in event.mimeData().urls():
            filename = str(url.encodedPath())[1:]
            print filename
            try:
                self._newSettings(filename)
            except:
                print "Unable to load settings file ", filename

    def _saveImageNumpy(self):
        self._control.saveFrameToNumpy("test_frame.npy")

    def _newSettings(self, filename):
        #self._control.saved_settings = [filename] + self._control.saved_settings
        self.settings = self._control.newSettings(filename)
        if len(self.settings) > self._max_saved:
            self.settings.pop()
        for i,p in enumerate(self.settings):
            filename = p["settings_filename"]
            radiobutton = getattr(self._ui, "settings"+str(i+1)+"_radioButton")
            radiobutton.setText(filename.split('/')[-1][:-5])
            radiobutton.show()

        self._ui.settings1_radioButton.click()

    def toggleSettings(self):
        for i,p in enumerate(self.settings):
            radioButton = getattr(self._ui, "settings"+str(i+1)+"_radioButton")
            if radioButton.isChecked():
                self._control.updateSettings(p)
        time.sleep(0.05)
        self._getDisplayTimings()
        self._writeROIValues()
        
        
    def _getDisplayTimings(self):
        self._exp = self._control.getExposureTime()
        self._ui.labelExposure.setText("Actual: %f ms" % (self._exp*1000))

    def _setTimings(self):
        exp_time = int(self._ui.lineEdit_ExposureTime.text())
        #delay_time = int(self._ui.lineEdit_DelayTime.text())
        exp_units = str(self._ui.comboBox_expunits.currentText())
        #delay_units = str(self._ui.comboBox_delayunits.currentText())

        if exp_units not in ['s','ms','us','ns']:
            exp_units = 's'

        '''
        if delay_units not in ['s','ms','us','ns']:
            delay_units = 's'
        '''

        if exp_units == 'ns':
            exp_time *= 1e-9
        elif exp_units == 'us':
            exp_time *= 1e-6
        elif exp_units == 'ms':
            exp_time *= 1e-3

        self._control.setExposureTime(exp_time)
        self._getDisplayTimings()

    def _resetPreview(self):
        self._previewMode = False
        self._control.stopCapture()
        self._control.beginPreview()
        self._previewMode = True

    def _setBinning(self):
        if self._ui.radioButton_2x2Binning.isChecked():
            self._control._set2x2Binning()
        else:
            self._control._setNoBinning()
        self._writeROIValues()

    def _writeROIValues(self):
        x0,y0 = self._control._getX0Y0()
        self._xres,self._yres = self._control._getResolution()
        self._ui.x0_label.setText("%i" % x0)
        self._ui.xres_label.setText("%i" % self._xres)
        self._ui.y0_label.setText("%i" % y0)
        self._ui.yres_label.setText("%i" % self._yres)
        self._ui.labelResolution.setText("%i, %i" % (self._xres,self._yres))
        x1 = x0+self._xres
        y1 = y0+self._yres
        self._ROIDisplayHelper(x0,x1,y0,y1)
        

    def _setROI(self):
        self._control.stopCapture()
        x0 = int(self._ui.x0_lineEdit.text())
        y0 = int(self._ui.y0_lineEdit.text())
        xres = int(self._ui.xres_lineEdit.text())
        yres = int(self._ui.yres_lineEdit.text())
        if x0>=0 and x0<2048 and y0>=0 and y0<2048:
            self._control._setROI(x0,y0,xres,yres)
            time.sleep(0.1)
            self._resetPreview()
        self._writeROIValues()

    def _setROItoCommon(self, setting):
        self._control.stopCapture()
        self._control.setROI_commonsettings(setting)
        self._writeROIValues()
        time.sleep(0.1)
        self._resetPreview()

    def _setPlot(self, dimension):
        if dimension == 'w':
            self.plotWidth = int(self._ui.lineEdit_plotWidth.text())
        elif dimension == 'l':
            self.plotLength = int(self._ui.lineEdit_plotLength.text())
        elif dimension == 'x':
            self.xForPlot = int(self._ui.lineEdit_xForPlot.text())
        elif dimension == 'y':
            self.yForPlot = int(self._ui.lineEdit_yForPlot.text())
        

    def _ROIDisplayHelper(self, x0,x1,y0,y1):
        roi_display = np.zeros((320,320),dtype=np.uint8)
        x0_scaled = int(x0/6.4)
        x1_scaled = int(x1/6.4)
        y0_scaled = int(y0/6.4)
        y1_scaled = int(y1/6.4)
        print x0_scaled, x1_scaled
        print y0_scaled, y1_scaled
        roi_display[y0_scaled:y1_scaled+1, x0_scaled:x1_scaled+1] = 255/2
        qt_roiimage = qext.numpy_to_qimage8(roi_display, 0,255,self._cmap)
        roi_pixmap = QtGui.QPixmap.fromImage(qt_roiimage)
        roi_pixmap = roi_pixmap.scaled(320/3, 320/3, QtCore.Qt.KeepAspectRatio)
        self._ui.labelROIDisplay.setPixmap(roi_pixmap)

            

    def _update(self):
        if self._previewMode:
            np_image = self._control.getImageForPreview()
        elif self._recordingOn:
            np_image = self._control.imageForPreviewWhileRecording()
        else:
            np_image = None
        if np_image is not None:
            self._plotCrossSection(np_image)
            np_min = np_image[2:-2,2:-2].min()
            np_max = np_image[2:-2,2:-2].max()
            if self._autoscale:
                self.vmin = np_min
                self.vmax = np_max
                self._ui.horizontalSlider_min.setValue(int(np_min))
                self._ui.horizontalSlider_max.setValue(int(np_max))
            if self._autoscale90:
                self.vmin = np_min
                self.vmax = 0.75*np_max
                self._ui.horizontalSlider_min.setValue(int(np_min))
                self._ui.horizontalSlider_max.setValue(int(np_max*0.75))
            if self.vmin==self.vmax or self.vmin>self.vmax:
                self.vmax = 1+self.vmin
            qt_image = qext.numpy_to_qimage8(np_image, self.vmin, self.vmax, self._cmap)
            self._pixmap = QtGui.QPixmap.fromImage(qt_image)
            #xdim = np.minimum(512, self._control._props["dimensions"][0])
            #ydim = np.minimum(512, self._control._props["dimensions"][1])
            #self._pixmap = self._pixmap.scaled(xdim, ydim, QtCore.Qt.KeepAspectRatio)
            self._pixmap = self._pixmap.scaled(512, 512, QtCore.Qt.KeepAspectRatio)
            #if self._ui.checkBox_showTarget.isChecked():
            #    self._drawTarget()
            self._ui.labelDisplay.update()
            self._ui.labelMin.setText(str(np_min))
            self._ui.labelMax.setText(str(np_max))
            self._ui.labelMean.setText(str(int(np_image.mean())))
            self._ui.labelMedian.setText(str(int(np.median(np_image))))
            #self._ui.labelDisplay.setPixmap(self._pixmap)
            self._ui.label_xy.setText("x,y: %i, %i; value: %i" % (self.x, self.y, np_image[self.y,self.x]))
            small_roi = np_image[self.y-5:self.y+5,self.x-5:self.x+5]
            if small_roi.shape[0]>0:
                if small_roi.shape[1]>0:
                    self._ui.label_xyArea.setText("%i" % np.max(np_image[self.y-5:self.y+5,self.x-5:self.x+5]))

    def _labelDisplay_paintEvent(self, event):
        qp = QtGui.QPainter(self._ui.labelDisplay)
        if self._pixmap != None:
            qp.drawPixmap(0,0,self._pixmap)
        qp.setPen(QtCore.Qt.red)
        if self._ui.checkBox_showTarget.isChecked():
            qp.drawLine(251,251,253,253)
            qp.drawLine(262,262,260,260)
            qp.drawLine(253,260,251,262)
            qp.drawLine(262,251,260,253)

    def _plotCrossSection(self, image):
        if self.plotWidth>0 and self.plotLength>0:
            hor = self._ui.radioButton_plotHor.isChecked()
            area = image[self.xForPlot:self.xForPlot+self.plotLength,
                         self.yForPlot:self.yForPlot+self.plotWidth]
            if hor:
                crossSection = area.sum(axis=1)
            else:
                crossSection = area.sum(axis=0)
            xaxis = np.arange(0,len(crossSection))
            curve1 = Qwt5.QwtPlotCurve('')
            pen = QtGui.QPen(QtCore.Qt.black)
            pen.setStyle(QtCore.Qt.SolidLine)
            curve1.setPen(pen)
            curve1.attach(self._ui.qwtPlot)
            curve1.setData(xaxis,crossSection)
            self._ui.qwtPlot.replot()
            curve1.detach()
        #if self._ui.checkBox_sharpness.isChecked():
            
            

    def _setFilename(self):
        filenum = "%.4i" % self._ui.spinBox_fileNameIndex.value()
        self.filename = self._ui.lineEdit_fileName.text() + filenum
        if os.path.exists(self.filename+".dax"):
            self._ui.label_fn.setStyleSheet("QLabel { color: red}")
            self._ui.spinBox_fileNameIndex.setStyleSheet("QLabel { color: red}")
            self._ui.spinBox_fileNameIndex.setStyleSheet("QLabel { border-color: red}")
        else:
            self._ui.label_fn.setStyleSheet("QLabel { color: black}")
            self._ui.spinBox_fileNameIndex.setStyleSheet("QLabel { color: black}")
            self._ui.spinBox_fileNameIndex.setStyleSheet("QLabel { border-color: black}")

    def _setRecordMode(self, mode):
        mode = self._ui.comboBox_recordMode.currentIndex()
        self._ui.label_fixedLength.setEnabled(bool(mode))
        self._ui.lineEdit_fixedLength.setEnabled(bool(mode))

    def _record(self):
        if not self._recordingOn:
            if self._ui.lineEdit_fixedLength.isEnabled():
                numToRecord = int(self._ui.lineEdit_fixedLength.text())
                self._control._setRecordMode(1, numToRecord)
            else:
                self._control._setRecordMode(0, 100)
            self._ui.pushButton_Record.setText('STOP')
            self._ui.pushButton_Record.setStyleSheet("QPushButton { color: red }")
            self._previewMode = False
            self._ui.comboBox_recordMode.setEnabled(False)
            self._control.beginRecording(self.filename)
            if self._ui.checkBox_autoUpdateTime.isChecked():
                uptime = (self._exp * 1000) + 10
                if (1000*self._exp)>180:
                    uptime = 400
                if uptime<500 and uptime>20:
                    self._toDiskUpdateTime = uptime
                    self._ui.lineEdit_toDiskTime.setText(str(uptime))
            self._updateToDisk.start(self._toDiskUpdateTime)
            self._recordingOn=True
        else:
            self._recordStop()

    def _recordStop(self):
        self._ui.pushButton_Record.setText('Record')
        self._ui.pushButton_Record.setStyleSheet("QPushButton { color: black }")
        self._ui.spinBox_fileNameIndex.setValue(self._ui.spinBox_fileNameIndex.value()+1)
        self._ui.comboBox_recordMode.setEnabled(True)
        self._updateToDisk.stop()
        self._control.endRecording()
        self._recordingOn=False
        self._control.stopCapture()
        self._control.beginPreview()
        self._previewMode = True
            

    def _fillDAX(self):
        nIms, nImsRec, done = self._control.recordToDAX()
        self._ui.labelImagesRecorded.setText(str(nIms) + ", " + str(nImsRec))
        if done:
            self._recordStop()

    def _beginSnaps(self):
        self._ui.pushButton_Record.setEnabled(False)
        self._control.beginSnaps(self.filename)
        self._recordingOn=True

    def _takeSnap(self):
        imNum = self._control.takeSnap()
        self._ui.labelImagesRecorded.setText(str(imNum))

    def _doneSnaps(self):
        self._control.endRecording()
        self._recordingOn=False
        self._ui.pushButton_Record.setEnabled(True)
        self._ui.spinBox_fileNameIndex.setValue(self._ui.spinBox_fileNameIndex.value()+1)

    def stopCapture(self):
        self._control.stopCapture()

    def setAutoscale(self, state):
        self._autoscale = state
        if state:
            self._ui.checkBox_autoscale90.setChecked(False)

    def setAutoscale90(self,state):
        self._autoscale90 = state
        if state:
            self._ui.checkBoxAutoscale.setChecked(False)


    def shutDown(self):
        self._updater.stop()
        self._updateToDisk.stop()
