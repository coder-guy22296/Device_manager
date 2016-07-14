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

class popUpOptions(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        vbox = QtGui.QVBoxLayout()
        vbox.addStretch(1)
        self.setWindowTitle('Additional Recording Options')
        #self.checkBox_killLasers = QtGui.QCheckBox("After recording, kill lasers?", self)
        #self.checkBox_killLasers.stateChanged.connect(self.updateKillLasers)
        self.checkBox_stopFL = QtGui.QCheckBox("After recording, stop focus lock?", self)
        self.checkBox_stopFL.setChecked(True)
        self.checkBox_stopFL.stateChanged.connect(self.updateStopFL)
        self.checkBox_enableShuttersOnStart = QtGui.QCheckBox("On record, ready shutters?", self)
        self.checkBox_enableShuttersOnStart.setChecked(True)
        self.checkBox_enableShuttersOnStart.stateChanged.connect(self.updateEnableShutters)
        #vbox.addWidget(self.checkBox_killLasers)
        vbox.addWidget(self.checkBox_stopFL)
        vbox.addWidget(self.checkBox_enableShuttersOnStart)
        

        self.killLasers = False
        self.stopFL = True
        self.readyShuttersOnRecord = True

        self.setGeometry(20,20,400,200)
        self.setLayout(vbox)

    def updateKillLasers(self):
        self.killLasers = self.checkBox_killLasers.isChecked()

    def updateStopFL(self):
        self.stopFL = self.checkBox_stopFL.isChecked()

    def updateEnableShutters(self):
        self.readyShuttersOnRecord = self.checkBox_enableShuttersOnStart.isChecked()
        
class spotFinder(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        '''
        self.view = pg.GraphicsView()
        self.graph = pg.PlotItem()
        self.view.setCentralWidget(self.graph)
        self.setCentralWidget(self.view)

        testx = np.random.rand(10)
        testy = np.random.rand(10)

        self.graph.plot(testx,testy,pen=None,symbol='o')
        '''

        self.qwtPlot = Qwt5.QwtPlot(self)
        self.setCentralWidget(self.qwtPlot)
        #self.qwtPlot.setGeometry(QtCore.QRect(260, 400, 581, 200))
        self.scatter = Qwt5.QwtPlotCurve('')
        self.scatter.attach(self.qwtPlot)
        self.scatter.setPen(QtGui.QPen(QtCore.Qt.NoPen))
        scatter_symbol = Qwt5.QwtSymbol(Qwt5.QwtSymbol.Ellipse,
                                        QtGui.QBrush(QtCore.Qt.red),
                                        QtGui.QPen(QtCore.Qt.red),
                                        QtCore.QSize(7,7))
        self.scatter.setSymbol(scatter_symbol)

        self.scatterAll = Qwt5.QwtPlotCurve('')
        self.scatterAll.attach(self.qwtPlot)
        self.scatterAll.setPen(QtGui.QPen(QtCore.Qt.NoPen))
        scatter_symbol = Qwt5.QwtSymbol(Qwt5.QwtSymbol.Ellipse,
                                        QtGui.QBrush(QtCore.Qt.black),
                                        QtGui.QPen(QtCore.Qt.black),
                                        QtCore.QSize(3,3))
        self.scatterAll.setSymbol(scatter_symbol)
                                           
                                           
        
        self.setGeometry(40,40,512,512)

        self.cell_size = 32
        self.threshold = 200

        self.fof = fof.MedFastObjectFinder(self.cell_size, self.threshold)
        #self.setLayout(vbox)

        self.statusBar = QtGui.QStatusBar(self)
        self.setStatusBar(self.statusBar)

        self.allXs = np.array([])
        self.allYs = np.array([])

    def setFoundSpots(self,n):
        self.statusBar.showMessage("Found spots in last frame: %i" % n)

    def setThreshold(self, threshold):
        self.threshold = threshold
        self.fof.changeThreshold(self.threshold)

    def setCellSize(self, size):
        self.cell_size = size
        self.fof.changeCellSize(self.cell_size)

    def updatePlot(self,x,y):
        self.allXs = np.hstack((self.allXs,x))
        self.allYs = np.hstack((self.allYs,y))
        #self.graph.plot(x,y,pen=None,symbol='o')
        self.scatter.setData(x,-1*y)
        self.scatterAll.setData(self.allXs,-1*self.allYs)
        self.qwtPlot.replot()

    def resetSpots(self, thresh, size):
        self.setThreshold(thresh)
        self.setCellSize(size)
        self.allXs = np.array([])
        self.allYs = np.array([])
        


class UI(inLib.ModuleUI):

    def __init__(self, control, ui_control):
        inLib.ModuleUI.__init__(self, control, ui_control, 'modules.STORMimaging.STORMimaging_design')

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

        self.additionalOptions = popUpOptions()
        self._ui.pushButton_moreOptions.clicked.connect(self.additionalOptions.show)

        self.spotFinder = spotFinder()
        self.spotFinder.show()
        self._ui.spotFinder_checkBox.clicked.connect(self.setSpotFinder)
        self.useSpotCounter=False
        self._ui.pushButton_resetSpots.clicked.connect(self.resetSpots)
        self._ui.lineEdit_threshold.setText("200")
        self._ui.lineEdit_cellsize.setText("32")

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

    def setSpotFinder(self):
        self.useSpotCounter = self._ui.spotFinder_checkBox.isChecked()

    def resetSpots(self):
        threshold = int(self._ui.lineEdit_threshold.text())
        cellsize = int(self._ui.lineEdit_cellsize.text())
        self.spotFinder.resetSpots(threshold, cellsize)
        

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
            if self.useSpotCounter:
                np_image_new = np_image.copy()
                np_image_new[0:5,:] = np_image.min()
                np_image_new[:,0:5] = np_image.min()
                np_image_new[:,-5:] = np_image.min()
                np_image_new[-5:,:] = np_image.min()
                frame = np_image_new.astype(np.int16).tostring()
                x,y,n = self.spotFinder.fof.findObjects(frame,np_image.shape[0],np_image.shape[1])
                xs = np.fromstring(x,np.dtype(np.float32))
                ys = np.fromstring(y,np.dtype(np.float32))
                self.spotFinder.setFoundSpots(n)
                self.spotFinder.updatePlot(xs,ys)
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
                if self.additionalOptions.readyShuttersOnRecord:
                    self._ui_control.illumination._ui.pushButton_readyShutters.click()
                    time.sleep(0.05)
                print "Running shutters..."
                self._ui_control.illumination._ui.pushButton_startShutters.click()
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
        if self.additionalOptions.stopFL:
            if self._ui_control.imageBasedFocusLock._ui.fl_pushButton.text() == 'Focus Lock Stop':
                self._ui_control.imageBasedFocusLock._ui.fl_pushButton.click()
        self._control._control.camera.endRecording()
        self._recordingOn=False
        self._control._control.camera.stopCapture()
        self._control._control.camera.beginPreview()
        self._previewMode = True
        #if self.additionalOptions.killLasers:
        self._control.offLasers()
        for onOffBoxes in self._ui_control.illumination.onoff:
            onOffBoxes.setCheckState(0)
        self._ui_control.illumination._ui.On561.setCheckState(0)
        self._ui_control.illumination._ui.OnLED.setCheckState(0)
        self._ui_control.illumination.stopShutters()
        
        #time.sleep(0.1)
        #self._control._control.shutters.stopDO()
            
            #print self._ui_control.illumination
            

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

    def showNewWindow(self):
        self.additionalOptions.setGeometry(QtCore.QRect(20,20,400,100))
        self.additionalOptions.show()


    def shutDown(self):
        self._updater.stop()
        self._updateToDisk.stop()
