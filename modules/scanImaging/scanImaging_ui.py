#!/usr/bin/python

import inLib
from PyQt4 import QtCore, QtGui, Qwt5
from Utilities import QExtensions as qext
from functools import partial
from libs.objectFinder import fastObjectFinder as fof
import time
import numpy as np
import os
import pyqtgraph as pg

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
        inLib.ModuleUI.__init__(self, control, ui_control, 'modules.scanImaging.scanImaging_design')

        self._control._control.camera.beginPreview()

        self._autoscale = True
        self._autoscale90 = False
        self._pixmap = None
        self._cmap = None
        self.vmin = 0.0
        self.vmax = 65535

        self.filename = 'movie0001'

        self._previewMode = True

        self.nIms = 0
        self.modNum = 0
        self.zerothFrame = 0
        
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
        self._ui.comboBox_filename.currentIndexChanged.connect(self._setFilename)
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
        self._ui.lineEdit_toDiskTime.returnPressed.connect(self._newUpdateTime)
        self._ui.pushButton_roi128128.clicked.connect(partial(self._setROItoCommon, "128x128"))
        self._ui.pushButton_roi256256.clicked.connect(partial(self._setROItoCommon, "256x256"))
        self._ui.pushButton_roi512256.clicked.connect(partial(self._setROItoCommon, "512x256"))
        self._ui.pushButton_roi512512.clicked.connect(partial(self._setROItoCommon, "512x512"))
        self._ui.pushButton_spiralSnaps.clicked.connect(self._spiralSnaps)
        self._ui.loadSettings_pushButton.clicked.connect(self.loadSettingsFile)
        self._ui.setModNum_pushButton.clicked.connect(self.getModNum)
        self._ui.loadPositions_pushButton.clicked.connect(self.loadPositions)
        self._ui.savePositions_pushButton.clicked.connect(self.savePositions)
        

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

        #Combobox for filename suffix
        self._ui.comboBox_filename.addItem("")
        self._ui.comboBox_filename.addItem("_647")
        self._ui.comboBox_filename.addItem("_561")
        self._ui.comboBox_filename.addItem("_532")
        self._ui.comboBox_filename.addItem("_488")
        self._ui.comboBox_filename.addItem("_conv")
        self._ui.comboBox_filename.addItem("_LED")
        self._ui.comboBox_filename.setCurrentIndex(0)

        #Setting text:
        xdim,ydim = self._control._control.camera._getResolution()
        self._ui.labelResolution.setText("%i, %i" % (xdim,ydim))
        self._control._control.camera._props['dimensions'] = xdim,ydim
        self._xres, self._yres = xdim,ydim
        self._ui.lineEdit_fileName.setText("movie")

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

        #self.settings = [self._control._control.camera.loadedSettings]
        #print self.settings[0]
        self.settings = self._control.newSettings(None)
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


        self._ui.loc1_radioButton.hide()
        self._ui.loc2_radioButton.hide()
        self._ui.loc3_radioButton.hide()
        self._ui.loc4_radioButton.hide()
        self._ui.loc5_radioButton.hide()
        self._ui.loc6_radioButton.hide()
        self._ui.loc7_radioButton.hide()
        self._ui.loc8_radioButton.hide()
        self._ui.remember_pushButton.clicked.connect(self.rememberLocation)
        self._ui.clearLocations_pushButton.clicked.connect(self.clearLocations)
        self._ui.loc1_radioButton.toggled.connect(partial(self.goToLocation,1))
        self._ui.loc2_radioButton.toggled.connect(partial(self.goToLocation,2))
        self._ui.loc3_radioButton.toggled.connect(partial(self.goToLocation,3))
        self._ui.loc4_radioButton.toggled.connect(partial(self.goToLocation,4))
        self._ui.loc5_radioButton.toggled.connect(partial(self.goToLocation,5))
        self._ui.loc6_radioButton.toggled.connect(partial(self.goToLocation,6))
        self._ui.loc7_radioButton.toggled.connect(partial(self.goToLocation,7))
        self._ui.loc8_radioButton.toggled.connect(partial(self.goToLocation,8))
        

        self._ui.pushButton_overlay.clicked.connect(self.useOverlay)


        ##Handling file drops
        
        self._ui.tabWidget.__class__.dragEnterEvent = self._dragEnterEvent
        self._ui.tabWidget.__class__.dropEvent = self._dropEvent
        self._ui.tabWidget.setAcceptDrops(True)
        #self._ui.tabWidget.setDropEnabled(True)
        #self._ui.tabWidget.setDragEnabled(True)
        
        self._ui.labelDisplay.setGeometry(10,10,512,512)

        self._ui.labelDisplay.__class__.dragEnterEvent = self._dragEnterEvent
        self._ui.labelDisplay.__class__.dropEvent = self._dropEvent

        #For testing purposes. Hide if not in testing mode
        self._ui.groupBox_testing.hide()

        # A timer to update the image:
        self._updater = QtCore.QTimer()
        self._updater.timeout.connect(self._update)
        self._updater.start(50)

        self._updaterForStage = QtCore.QTimer()
        self._updaterForStage.timeout.connect(self._updateStageStatus)
        self._updaterForStage.start(300)

        self._updateToDisk = QtCore.QTimer()
        self._updateToDisk.timeout.connect(self._fillDAX)
        self._toDiskUpdateTime = 200

        #Updater for spiral snaps:
        self._spiralUpdate = QtCore.QTimer()
        self._spiralUpdate.timeout.connect(self._spiral)
        self._spiralNum = 0

        self._recordingOn = False

        self._ui.spinBox_fileNameIndex.setValue(1)

        self._ui.horizontalSlider_min.setRange(0,65535)
        self._ui.horizontalSlider_max.setRange(0,65535)
        self._ui.horizontalSlider_min.valueChanged.connect(self._setDispMin)
        self._ui.horizontalSlider_max.valueChanged.connect(self._setDispMax)

        self._ui.defectCorrect_checkBox.setChecked(True)
        self._ui.defectCorrect_checkBox.stateChanged.connect(self.setDefectCorrectMode)

        self._ui.exposureOutputTrigger_pushButton.clicked.connect(self.expOutTrig)
        self._ui.progOutputTrigger_pushButton.clicked.connect(self.progOutTrig)

        self.killLasers=False
        self.overlay2Color = False

        self.posRadioButtons = [self._ui.pos1_radioButton,
                                self._ui.pos2_radioButton,
                                self._ui.pos3_radioButton,
                                self._ui.pos4_radioButton,
                                self._ui.pos5_radioButton,
                                self._ui.pos6_radioButton,
                                self._ui.pos7_radioButton,
                                self._ui.pos8_radioButton,
                                self._ui.pos9_radioButton,
                                self._ui.pos10_radioButton]
        i=0
        for button in self.posRadioButtons:
            button.hide()
            button.clicked.connect(partial(self._goToStoredLoc, i))
            i=i+1

        self._ui.storePos_pushButton.clicked.connect(self._storePos)
        self._ui.clearPos_pushButton.clicked.connect(self._clearPos)

        self._ui.startBackAndForth_pushButton.clicked.connect(self.startBackAndForth)
        self._ui.stopBackAndForth_pushButton.clicked.connect(self.stopBackAndForth)
        self._ui.scan3pos_pushButton.clicked.connect(self.start3pos)
        self._ui.scan1to2wait_lineEdit.returnPressed.connect(self._setWaitTime)
        self._ui.getStageStatus_pushButton.clicked.connect(self.getStageStatus)

        self._ui.blockLaser_pushButton.clicked.connect(self.blockLaser)
        self._ui.unblockLaser_pushButton.clicked.connect(self.unblockLaser)

        self._ui.startPump_pushButton.clicked.connect(self.startPump)
        self._ui.stopPump_pushButton.clicked.connect(self.stopPump)
        self._ui.pumpFlowRate_lineEdit.returnPressed.connect(self.pumpFlowRate)

        self._backAndForthTimer = QtCore.QTimer()
        self._backAndForthTimer.timeout.connect(self.backAndForth)
        self._backAndForthTimer3 = QtCore.QTimer()
        self._backAndForthTimer3.timeout.connect(self.backAndForth3)
        self.current = 0

        self._ui.startScan_pushButton.clicked.connect(self.startScanRecord)
        self._ui.stopScan_pushButton.clicked.connect(self.stopScanRecord)
        self._scanRecord = None

        self.newSpeed = 0.05


    def startBackAndForth(self):
        self.newSpeed = float(self._ui.scan1to2speed_lineEdit.text())
        wtime = self._control.returnWaitTime()
        print "Wait time: ", wtime
        self._backAndForthTimer.start(wtime*1000)

    def start3pos(self):
        self.newSpeed = float(self._ui.scan1to2speed_lineEdit.text())
        wtime = self._control.returnWaitTime()
        self._backAndForthTimer3.start(wtime*1000)

    def stopBackAndForth(self):
        self._backAndForthTimer.stop()
        self._backAndForthTimer3.stop()

    def blockLaser(self):
        self._control.writePositionToArduino(20)

    def unblockLaser(self):
        self._control.writePositionToArduino(100)

    def _whatPumpNum(self):
        if self._ui.pump1_radioButton.isChecked():
            num = 1
        elif self._ui.pump0_radioButton.isChecked():
            num = 0
        else:
            num = 2
        return num

    def startPump(self):
        num = self._whatPumpNum()
        self._control.runPumpForward(num)

    def stopPump(self):
        num = self._whatPumpNum()
        self._control.stopPump(num)

    def pumpFlowRate(self):
        num = self._whatPumpNum()
        rate = int(self._ui.pumpFlowRate_lineEdit.text())
        self._control.setPumpForwardRate(num,rate)
        

    def _setWaitTime(self):
        wtime = float(self._ui.scan1to2wait_lineEdit.text())
        print "Setting wait time to ", wtime
        self._control.setWaitTime(wtime)

    def _goToStoredLoc(self,num):
        self._control.goToLocations(num)

    def _storePos(self):
        locs = self._control.addLocation()
        for i in range(0,len(locs)):
            self.posRadioButtons[i].setText("Pos %i: %.2f,%.2f" % (i,locs[i][0]/10000.,locs[i][1]/10000.))
            self.posRadioButtons[i].setVisible(True)

    def _clearPos(self):
        self._control.clearLocations()
        for button in self.posRadioButtons:
            button.hide()
            button.setText("Pos: --,--")

    def savePositions(self):
        filename = QtGui.QFileDialog.getSaveFileName(None, 'Save positions',
                                                     '','*.txt')
        self._control.savePositions(str(filename))


    def loadPositions(self):
        filename = QtGui.QFileDialog.getOpenFileName(None, 'Load positions',
                                                     '', '*.txt')
        if filename:
            self._clearPos()
            locs = self._control.loadPositions(str(filename))
            for i in range(0,len(locs)):
                self.posRadioButtons[i].setText("Pos %i: %.2f,%.2f" % (i,locs[i][0]/10000.,locs[i][1]/10000.))
                self.posRadioButtons[i].setVisible(True)
            
        

    def getStageStatus(self):
        resp = self._control.getStageStatus()

    def backAndForth(self):
        self._control.setSpeed('xy', self.newSpeed)
        self._control.goToLocations(int(self.current))
        self.current = int(not self.current)
        if self.current:
            self.newSpeed = float(self._ui.scan1to2speed_lineEdit.text())
        else:
            self.newSpeed = float(self._ui.scan2to1speed_lineEdit.text())
        

    def backAndForth3(self):
        self._control.setSpeed('xy', self.newSpeed)
        self._control.goToLocations(int(self.current))
        if self.current:
            self.current = 2
            self.newSpeed = float(self._ui.scan2to3speed_lineEdit.text())
        elif self.current == 2:
            self.current = 0
            self.newSpeed = float(self._ui.scan2to1speed_lineEdit.text())
        elif self.current == 0:
            self.current = 1
            self.newSpeed = float(self._ui.scan1to2speed_lineEdit.text())

    def _printCamProps(self):
        self._control._control.camera.printProperties()

    def _newUpdateTime(self):
        uptime = int(self._ui.lineEdit_toDiskTime.text())
        if uptime>10 and uptime<10000:
            self._toDiskUpdateTime = uptime

    def expOutTrig(self):
        self._control._control.camera.setToExposureMode()

    def progOutTrig(self):
        self._control._control.camera.setToProgrammableMode()

    def setDefectCorrectMode(self):
        self._control._control.camera.stopCapture()
        mode = self._ui.defectCorrect_checkBox.isChecked()
        if mode:
            self._control._control.camera.setDefectCorrect(2)
        else:
            self._control._control.camera.setDefectCorrect(1)
        self._control._control.camera.beginPreview()
        

    def clearLocations(self):
        self._control.clearLocations()
        self._ui.loc1_radioButton.hide()
        self._ui.loc2_radioButton.hide()
        self._ui.loc3_radioButton.hide()
        self._ui.loc4_radioButton.hide()
        self._ui.loc5_radioButton.hide()
        self._ui.loc6_radioButton.hide()
        self._ui.loc7_radioButton.hide()
        self._ui.loc8_radioButton.hide()

    def rememberLocation(self):
        num = self._control.rememberLocation()
        radiobut = getattr(self._ui,"loc"+str(int(num))+"_radioButton")
        radiobut.show()

    def goToLocation(self, num):
        radiobut = getattr(self._ui,"loc"+str(int(num))+"_radioButton")
        if radiobut.isChecked():
            print "Moving to location ", str(num)
            self._control.goToLocation(int(num))

    def useOverlay(self):
        self.overlay2Color = not self.overlay2Color

    def _lowLight(self, state):
        if state:
            self._ui.horizontalSlider_min.setRange(0,10000)
            self._ui.horizontalSlider_max.setRange(0,10000)
        else:
            self._ui.horizontalSlider_min.setRange(0,65535)
            self._ui.horizontalSlider_max.setRange(0,65535)

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
        self._control._control.camera._setProp(propID,propVal)

    def _getProp(self):
        propID = int(self._ui.lineEdit_propID.text())
        propVal = self._control._control.camera._getProp(propID)
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
        self._control._control.camera.saveFrameToNumpy("test_frame.npy")

    def _newSettings(self, filename):
        #self._control.saved_settings = [filename] + self._control.saved_settings
        self.settings = self._control.newSettings(filename)
        if len(self.settings) > self._max_saved:
            self.settings.pop()
        for i,p in enumerate(self.settings):
            if p.has_key("settings_filename"):
                filename = p["settings_filename"]
                buttonText = filename.split('/')[-1][:-5]
            else:
                buttonText = "unknown"
            radiobutton = getattr(self._ui, "settings"+str(i+1)+"_radioButton")
            radiobutton.setText(buttonText)
            radiobutton.show()

        self._ui.settings1_radioButton.click()

    def loadSettingsFile(self):
        settingsFile = str(QtGui.QFileDialog.getOpenFileName(self._window,
                                                             'Open File',
                                                             'D:\\'))
        if settingsFile == '' or settingsFile is None:
            return
        else:
            self._newSettings(settingsFile)

    def toggleSettings(self):
        for i,p in enumerate(self.settings):
            radioButton = getattr(self._ui, "settings"+str(i+1)+"_radioButton")
            if radioButton.isChecked():
                self._control.updateSettings(p)
        time.sleep(0.05)
        self._getDisplayTimings()
        self._writeROIValues()
        
        
    def _getDisplayTimings(self):
        #self._exp = self._control._control.camera._getExposureTime()
        self._exp, fps = self._control.getExposureTimeFrameRate()
        self._ui.labelExposure.setText("Actual: %f ms" % (self._exp*1000))
        self._ui.labelFPS.setText("FPS: %.2f" % fps)

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
        self._control._control.camera.stopCapture()
        self._control._control.camera.beginPreview()
        self._previewMode = True

    def _setBinning(self):
        if self._ui.radioButton_2x2Binning.isChecked():
            self._control._control.camera._set2x2Binning()
        else:
            self._control._control.camera._setNoBinning()
        self._writeROIValues()

    def _writeROIValues(self):
        #x0,y0 = self._control._control.camera._getX0Y0()
        #self._xres,self._yres = self._control._control.camera._getResolution()
        x0,y0,self._xres,self._yres = self._control.getROI()
        self._ui.x0_label.setText("%i" % x0)
        self._ui.xres_label.setText("%i" % self._xres)
        self._ui.y0_label.setText("%i" % y0)
        self._ui.yres_label.setText("%i" % self._yres)
        self._ui.labelResolution.setText("%i, %i" % (self._xres,self._yres))
        x1 = x0+self._xres
        y1 = y0+self._yres
        self._ROIDisplayHelper(x0,x1,y0,y1)
        

    def _setROI(self):
        self._control._control.camera.stopCapture()
        x0 = int(self._ui.x0_lineEdit.text())
        y0 = int(self._ui.y0_lineEdit.text())
        xres = int(self._ui.xres_lineEdit.text())
        yres = int(self._ui.yres_lineEdit.text())
        if x0>=0 and x0<2048 and y0>=0 and y0<2048:
            self._control._control.camera._setROI(x0,y0,xres,yres)
            time.sleep(0.1)
            self._resetPreview()
        self._writeROIValues()

    def _setROItoCommon(self, setting):
        self._control._control.camera.stopCapture()
        self._control._control.camera.setROI_commonsettings(setting)
        self._writeROIValues()
        time.sleep(0.1)
        self._resetPreview()
        

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

            
    def getModNum(self):
        self.modNum = int(self._ui.modNum_spinBox.value())
        self.zerothFrame = int(self._ui.zerothFrame_spinBox.value())

    def _updateStageStatus(self):
        x,y,z = self._control.getXYZPosition()
        status = self._control.getStageStatus()
        #print status
        if self._control.getStageStatus():
            self._ui.labelStage.setStyleSheet("QLabel { color: black}")
        else:
            self._ui.labelStage.setStyleSheet("QLabel { color: red}")
        self._ui.labelStage.setText("%.3f, %.3f" % (x/10000., y/10000.))

    def _update(self):
        #xstage,ystage,zstage = self._control.getXYZPosition()
        #self._ui.labelStage.setText("%.1f, %.1f" % (xstage,ystage))
        if self._previewMode:
            np_image = self._control._control.camera.getImageForPreview()
        elif self._recordingOn:
            np_image = self._control._control.camera.imageForPreviewWhileRecording()
            if self.modNum>0:
                if np.mod(self.nIms-self.zerothFrame,self.modNum)>0:
                    np_image = None
        else:
            np_image = None
        if np_image is not None:
            if self.overlay2Color:
                np_image2 = np.zeros((256,256), dtype=np.uint16)
                np_image = np.rot90(np_image.reshape(256,512))
                np_image2 += (0.5*np_image[0:256,:].copy()) + (0.5*np_image[256:,:].copy())
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
            if self.overlay2Color:
                qt_image = qext.numpy_to_qimage8(np_image2, self.vmin, self.vmax, self._cmap)
            else:
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
            if np_image.shape[0]>self.x and np_image.shape[1]>self.y:
                self._ui.label_xy.setText("x,y: %i, %i; value: %.1f" % (self.x, self.y, np_image[self.y,self.x]))

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
        if self._ui.checkBox_showGrid.isChecked():
            qp.drawLine(64,0,64,511)
            qp.drawLine(128,0,128,511)
            qp.drawLine(196,0,196,511)
            qp.drawLine(256,0,256,511)
            qp.drawLine(320,0,320,511)
            qp.drawLine(384,0,384,511)
            qp.drawLine(448,0,448,511)
            qp.drawLine(0,64,511,64)
            qp.drawLine(0,128,511,128)
            qp.drawLine(0,196,511,196)
            qp.drawLine(0,256,511,256)
            qp.drawLine(0,320,511,320)
            qp.drawLine(0,384,511,384)
            qp.drawLine(0,448,511,448)

    def _setFilename(self):
        filenum = "%.4i" % self._ui.spinBox_fileNameIndex.value()
        suffix = self._ui.comboBox_filename.currentText()
        self.filename = self._ui.lineEdit_fileName.text() + '_' + filenum + suffix
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
            if self._ui.checkBox_runShutters.isChecked():
                self.stopCapture()
                time.sleep(0.05)
            if self._ui.lineEdit_fixedLength.isEnabled():
                numToRecord = int(self._ui.lineEdit_fixedLength.text())
                self._control._control.camera._setRecordMode(1, numToRecord)
            else:
                self._control._control.camera._setRecordMode(0, 100)
            self._ui.pushButton_Record.setText('STOP')
            self._ui.pushButton_Record.setStyleSheet("QPushButton { color: red }")
            self._previewMode = False
            self._ui.comboBox_recordMode.setEnabled(False)
            self._control._control.camera.beginRecording(self.filename)
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
        self._control._control.camera.endRecording()
        self._recordingOn=False
        self._control._control.camera.stopCapture()
        self._control._control.camera.beginPreview()
        self._previewMode = True
        #if self.additionalOptions.killLasers:
            

    def _fillDAX(self):
        nIms, nImsRec, done = self._control._control.camera.recordToDAX()
        self.nIms = nIms
        self._ui.labelImagesRecorded.setText(str(nIms) + ", " + str(nImsRec))
        if done:
            self._recordStop()

    def _beginSnaps(self):
        self._ui.pushButton_Record.setEnabled(False)
        self._control._control.camera.beginSnaps(self.filename)
        self._recordingOn=True

    def _takeSnap(self):
        imNum = self._control._control.camera.takeSnap()
        self._ui.labelImagesRecorded.setText(str(imNum))

    def _doneSnaps(self):
        self._control._control.camera.endRecording()
        self._recordingOn=False
        self._ui.pushButton_Record.setEnabled(True)
        self._ui.spinBox_fileNameIndex.setValue(self._ui.spinBox_fileNameIndex.value()+1)

    def _spiralSnaps(self):
        self._spiralNum = 0
        self._spiralUpdate.start(300)

    def _spiral(self):
        xbegin, ybegin, z = self._control.getXYZPosition()
        dx,dy = 1,1
        if self._spiralNum<100:
            if self._spiralNum>10:
                dy = 0
                dx = -1
            if self._spiralNum>30:
                dx,dy = -1,-1
            if self._spiralNum>50:
                dy = 0
                dx = 1
            if self._spiralNum>80:
                dx,dy = 1,1
            if self._spiralNum>95:
                dx,dy = 0,-1
            im = self._control.moveAndSnap(dx/10.0,dy/10.0)
            self._ui.labelImagesRecorded.setText(str(im))
            self._spiralNum += 1
            
        
    def stopCapture(self):
        self._control._control.camera.stopCapture()

    def setAutoscale(self, state):
        self._autoscale = state
        if state:
            self._ui.checkBox_autoscale90.setChecked(False)

    def setAutoscale90(self,state):
        self._autoscale90 = state
        if state:
            self._ui.checkBoxAutoscale.setChecked(False)

    def startScanRecord(self):
        startPos = self._ui.startScan_spinBox.value()
        stopPos = self._ui.stopScan_spinBox.value()
        restPos = self._ui.rest_spinBox.value()
        iterations = int(self._ui.scanIterations_lineEdit.text())
        scanSpeed = float(self._ui.scanSpeed_lineEdit.text())
        waitTime = float(self._ui.waitScanTime_lineEdit.text())
        start2 = self._ui.startScan2_spinBox.value()
        end2 = self._ui.stopScan2_spinBox.value()
        start3 = self._ui.startScan3_spinBox.value()
        end3 = self._ui.stopScan3_spinBox.value()
        start4 = self._ui.startScan4_spinBox.value()
        end4 = self._ui.stopScan4_spinBox.value()
        if self._ui.scan1_radioButton.isChecked():
            startPositions = [startPos]
            endPositions = [stopPos]
        elif self._ui.scan2_radioButton.isChecked():
            startPositions = [startPos, start2]
            endPositions = [stopPos, end2]
        elif self._ui.scan3_radioButton.isChecked():
            startPositions = [startPos, start2, start3]
            endPositions = [stopPos, end2, end3]
        elif self._ui.scan4_radioButton.isChecked():
            startPositions = [startPos, start2, start3, start4]
            endPositions = [stopPos, end2, end3, end4]
        else:
            startPositions = [startPos]
            endPositions = [stopPos]
        pumpAtRest = self._ui.onlyPumpOnRest_checkBox.isChecked()
        self._scanRecord = ScanAndRecord(self._control, self, startPositions,
                                         endPositions, restPos, scanSpeed, iterations,
                                         waitTime, pumpAtRest)
        self._window.connect(self._scanRecord, QtCore.SIGNAL('beginRecord'), self._record)
        self._window.connect(self._scanRecord, QtCore.SIGNAL('endRecord'), self._recordStop)
        print "self._scanRecord: ", self._scanRecord
        self._scanRecord.start()

    def stopScanRecord(self):
        self._scanRecord.stop()
        self._scanRecord = None
        


    def shutDown(self):
        self._updater.stop()
        self._updateToDisk.stop()

class ScanAndRecord(QtCore.QThread):
    def __init__(self, control, ui, startPos, end, rest_position,
                 scanSpeed, iterations, wTime, pumpAtRest):
        QtCore.QThread.__init__(self)

        self.control = control
        self.ui = ui
        self.startPos = startPos
        self.end = end
        self.rest_position = rest_position
        self.scanSpeed = scanSpeed
        self.iterations = iterations
        self.waitTime = wTime
        self.live = False
        self.fastScanSpeed = 1
        self.pumpAtRest = pumpAtRest

    def run(self):
        self.live=True
        for i in range(self.iterations):
            if self.live:
                for j in range(0,len(self.startPos)):
                    if self.pumpAtRest:
                        self.ui.stopPump()
                    self.control.goToLocations(self.startPos[j])
                    time.sleep(18)
                    self.control.setSpeed('xy', self.scanSpeed)
                    time.sleep(2)
                    self.control.goToLocations(self.end[j])
                    time.sleep(0.05)
                    self.emit(QtCore.SIGNAL('beginRecord'))
                    #self.ui._record()
                    while not self.control.getStageStatus():
                        time.sleep(0.5)
                    #self.ui._recordStop()
                    self.emit(QtCore.SIGNAL('endRecord'))
                    self.control.setSpeed('xy', self.fastScanSpeed)
                    self.control.goToLocations(self.rest_position)
                    if self.pumpAtRest:
                        self.ui.startPump()
                    time.sleep(self.waitTime*60)

    def stop(self):
        self.live = False
        print "Ending ScanAndRecord thread..."
