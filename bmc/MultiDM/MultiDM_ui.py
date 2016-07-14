#!/usr/bin/python


from PyQt4 import QtGui,QtCore
import inLib
from Utilities import QExtensions as qext
import numpy as np
from numpy.lib.scimath import sqrt as _msqrt
import copy
import time
import skimage

class UI(inLib.DeviceUI):
    
    def __init__(self, control):

        #path to design_path
        design_path = 'bmc.MultiDM.mirror_design'
        inLib.DeviceUI.__init__(self, control, design_path)

        self._ui.pushButton_load.clicked.connect(self.loadPattern)
        self._ui.pushButton_rot90.clicked.connect(self.patternRot90)
        self._ui.pushButton_fliplr.clicked.connect(self.patternFlipLR)
        self._ui.pushButton_flipud.clicked.connect(self.patternFlipUD)
        self._ui.pushButton_rotate.clicked.connect(self.patternRotate)
        self._ui.pushButton_getSegs.clicked.connect(self.getSegments)
        self._ui.pushButton_toMirror.clicked.connect(self.toMirror)
        self._ui.pushButton_toMirrorManyZerns.clicked.connect(self.toMirrorManyZerns)
        self._ui.pushButton_toMirrorManyZernRadii.clicked.connect(self.toMirrorManyZernRadii)
        self._ui.pushButton_toMirrorGroupVary.clicked.connect(self.toMirrorGroupVary)
        self._ui.pushButton_varyMult.clicked.connect(self.toMirrorVaryMult)
        self._ui.pushButton_reconfig.clicked.connect(self.reconfig)
        self._ui.pushButton_mult.clicked.connect(self.setMultiplier)
        self._ui.pushButton_premult.clicked.connect(self.setPreMultiplier)
        self._ui.pushButton_poke.clicked.connect(self.pokeSegment)
        self._ui.pushButton_clear.clicked.connect(self.clearPattern)
        self._ui.pushButton_refresh.clicked.connect(self.refreshPattern)
        self._ui.pushButton_pad.clicked.connect(self.padZeros)
        self._ui.pushButton_applyZern.clicked.connect(self.calcZernike)
        self._ui.pushButton_modulateZernike.clicked.connect(self.modZernike)
        self._ui.pushButton_createGroup.clicked.connect(self.createGroup)
        self._ui.pushButton_setToGroup.clicked.connect(self.setGroupVal)

        self._ui.pushButton_loadSegs.clicked.connect(self.loadSegs)

        self._ui.pushButton_loadvaryfile.clicked.connect(self.loadVaryFile)
        self.varyfilename = None
        self._ui.pushButton_varyfile.clicked.connect(self.varyfile)

        self._ui.lineEdit_loadMult.setText("10")
        self._ui.lineEdit_npixels.setText(str(self._control.pixels))
        self._ui.lineEdit_zernAmp.setText("0")
        self._ui.lineEdit_premult.setText(str(self._control.preMultiplier))

        self.pattern=None

        self._applyToMirrorThread = None
        self._applyManyZernsThread = None
        self._applyManyZernRadiiThread = None
        self._applyGroupOffsetsThread = None
        self._applyManyMultsToMirrorThread = None

    def loadPattern(self):
        filename = QtGui.QFileDialog.getOpenFileName(None,'Open pattern','','*.npy')
        m = float(self._ui.lineEdit_loadMult.text())
        pattern = self._control.loadPattern(filename, m)
        self._displayPhase(pattern)

    def loadVaryFile(self):
        self.varyfilename = QtGui.QFileDialog.getOpenFileName(None,'Open list of files file','','*.*')
        self._ui.label_filenameloaded.setText(self.varyfilename)

    def loadSegs(self):
        filename = QtGui.QFileDialog.getOpenFileName(None,'Open segments','','*.*')
        self._control.loadSegments(str(filename))
        segments = self._control.getSegments()
        self._displaySegments(segments)

    def clearPattern(self):
        self._control.clear()
        self._displayPhase(self._control.returnPattern())
        self.reportGeo()

    def refreshPattern(self):
        self._displayPhase(self._control.returnPattern())

    def patternRot90(self):
        pattern = self._control.patternRot90()
        self._displayPhase(pattern)

    def patternFlipLR(self):
        pattern = self._control.patternFlipLR()
        self._displayPhase(pattern)

    def patternFlipUD(self):
        pattern = self._control.patternFlipUD()
        self._displayPhase(pattern)

    def patternRotate(self):
        rot = float(self._ui.lineEdit_rotate.text())
        pattern = self._control.patternRotate(rot)
        self._displayPhase(pattern)

    def getSegments(self):
        self._control.findSegments()
        segments = self._control.getSegments()
        self._displaySegments(segments)

    def reconfig(self):
        cx = int(self._ui.lineEdit_cx.text())
        cy = int(self._ui.lineEdit_cy.text())
        npixels = int(self._ui.lineEdit_npixels.text())
        pattern = self._control.reconfigGeo(cx,cy,npixels)
        self._displayPhase(pattern)
        self.reportGeo()

    def setMultiplier(self):
        mult = float(self._ui.lineEdit_mult.text())
        self._control.setMultiplier(mult)

    def setPreMultiplier(self):
        mult = float(self._ui.lineEdit_premult.text())
        self._control.setPreMultiplier(mult)

    def toMirror(self):
        self._applyToMirrorThread = ApplyToMirror(self._control)
        self._applyToMirrorThread.start()
        #self._control.applyToMirror()

    def toMirrorManyZerns(self):
        '''
        runs when pushButton_toMirrorManyZerns is pushed

        Creates instance of class ApplyManyZernsToMirror which runs as
        a separate thread.
        '''
        mode = self._ui.spinBox_zernMode.value()
        mask = self._ui.checkBox_zernMask.isChecked()
        maxAmp = float(self._ui.lineEdit_maxZAmp.text())
        minAmp = float(self._ui.lineEdit_minZAmp.text())
        num = self._ui.spinBox_numZerns.value() #Zernike Mode
        wTime = int(self._ui.lineEdit_wTime.text())
        radius = int(self._ui.lineEdit_zernRad.text())
        '''
        if zernWithSharpness is checked, then advancing through zernike amplitudes
        will wait for AO module. Otherwise, will advance automatically
        '''
        withRunning = self._ui.checkBox_zernWithSharpness.isChecked()
        clearFirst = self._ui.checkBox_clearFirst.isChecked()
        self._applyManyZernsThread = ApplyManyZernsToMirror(self._control, mode, maxAmp,
                                                            minAmp, num, wTime, radius, mask,
                                                            clearFirst,
                                                            withRunning)
        self._applyManyZernsThread.start()

    def toMirrorManyZernRadii(self):
        mode = self._ui.spinBox_zernMode.value()
        amp = float(self._ui.lineEdit_zernAmp.text())
        mask = self._ui.checkBox_zernMask.isChecked()
        maxR = float(self._ui.lineEdit_maxZAmp.text())
        minR = float(self._ui.lineEdit_minZAmp.text())
        num = self._ui.spinBox_numZerns.value() #Zernike Mode
        wTime = int(self._ui.lineEdit_wTime.text())
        radius = int(self._ui.lineEdit_zernRad.text())
        '''
        if zernWithSharpness is checked, then advancing through zernike amplitudes
        will wait for AO module. Otherwise, will advance automatically
        '''
        withRunning = self._ui.checkBox_zernWithSharpness.isChecked()
        clearFirst = self._ui.checkBox_clearFirst.isChecked()
        self._applyManyZernRadiiThread = ApplyManyZernRadiiToMirror(self._control, mode, amp, maxR,
                                                                    minR, num, wTime, radius, mask,
                                                                    clearFirst,
                                                                    withRunning)
        self._applyManyZernRadiiThread.start()


    def toMirrorGroupVary(self):
        maxAmp = float(self._ui.lineEdit_maxZAmp.text())
        minAmp = float(self._ui.lineEdit_minZAmp.text())
        num = self._ui.spinBox_numZerns.value()
        wTime = int(self._ui.lineEdit_wTime.text())
        withRunning = self._ui.checkBox_zernWithSharpness.isChecked()
        clearFirst = self._ui.checkBox_clearFirst.isChecked()
        self._applyGroupOffsetsThread = ApplyManyOffsetsToGroup(self._control, maxAmp,
                                                                minAmp, num, wTime,
                                                                clearFirst,
                                                                withRunning)
        self._applyGroupOffsetsThread.start()

    def toMirrorVaryMult(self):
        maxMult = float(self._ui.lineEdit_maxZAmp.text())
        minMult = float(self._ui.lineEdit_minZAmp.text())
        num = self._ui.spinBox_numZerns.value()
        wTime = int(self._ui.lineEdit_wTime.text())
        withRunning = self._ui.checkBox_zernWithSharpness.isChecked()
        clearFirst = self._ui.checkBox_clearFirst.isChecked()

        self._applyManyMultsToMirrorThread = ApplyManyMultsToMirror(self._control, minMult,
                                                                    maxMult, num, wTime, withRunning)
        self._applyManyMultsToMirrorThread.start()

    def varyfile(self):
        wTime = int(self._ui.lineEdit_wTime.text())
        withRunning = self._ui.checkBox_zernWithSharpness.isChecked()
        if self.varyfilename is not None:
            self._applyVaryFileThread = ApplyManyFromFile(self._control,
                                                          self.varyfilename, wTime,
                                                          withRunning)
            self._applyVaryFileThread.start()

    def pokeSegment(self):
        segment = self._ui.spinBox_segment.value()
        toAdd = int(self._ui.lineEdit_pokeval.text())
        pokeAll = self._ui.checkBox_pokeAll.isChecked()
        self._control.pokeSegment(segment-1,toAdd,pokeAll=pokeAll)
        self._displayPhase(self._control.returnPattern())
        self._displaySegments(self._control.returnSegments())

    def padZeros(self):
        toPad = int(self._ui.lineEdit_pad.text())
        pattern = self._control.padZeros(toPad)
        self._displayPhase(pattern)
        self.reportGeo()

    def reportGeo(self):
        npixels, cx, cy = self._control.getGeoParams()
        self._ui.lineEdit_npixels.setText(str(npixels))
        self._ui.lineEdit_cx.setText(str(int(cx)))
        self._ui.lineEdit_cy.setText(str(int(cy)))

    def calcZernike(self):
        mode = self._ui.spinBox_zernMode.value()
        amplitude = float(self._ui.lineEdit_zernAmp.text())
        mask = self._ui.checkBox_zernMask.isChecked()
        radius = int(self._ui.lineEdit_zernRad.text())
        zern = self._control.calcZernike(mode, amplitude, radius=radius, useMask=mask)
        self._displayZern(zern)

    def modZernike(self):
        pattern = self._control.addZernike()
        self._displayPhase(pattern)

    def createGroup(self):
        groupStr = self._ui.lineEdit_group.text()
        group = np.array([int(s) for s in groupStr.split(',')]) - 1
        p = self._control.highlight_dummy_mirror_segs(group)
        self._displayDummyPattern(p)
        return group

    def setGroupVal(self):
        group = self.createGroup()
        toAdd = int(self._ui.lineEdit_groupVal.text())
        self._control.pokeGroup(group, toAdd)
        #for g in group:
        #    self._control.pokeSegment(g,toAdd,pokeAll=False)
        self._displayPhase(self._control.returnPattern())
        self._displaySegments(self._control.returnSegments())

    
    def _displayDummyPattern(self, pattern):
        if pattern is not None:
            self._ui.mplwidgetGrouped.figure.axes[0].matshow(pattern, cmap='RdBu')
            self._ui.mplwidgetGrouped.draw()

    def _displayZern(self, zernike):
        if zernike is not None:
            self._ui.mplwidgetZern.figure.axes[0].matshow(zernike, cmap='RdBu')
            self._ui.mplwidgetZern.draw()
        
    def _displayPhase(self, phase):
        if phase is not None:
            self._ui.mplwidgetPhase.figure.axes[0].matshow(phase, cmap='RdBu')
            self._ui.mplwidgetPhase.draw()

    def _displaySegments(self, segs):
        if segs is not None:
            self._ui.mplwidgetSegs.figure.axes[0].matshow(segs, cmap='RdBu')
            self._ui.mplwidgetSegs.draw()
        flatsegs = segs.flatten()
        trueSegs = np.zeros((140))
        trueSegs[0:10] = flatsegs[1:11]
        trueSegs[10:130] = flatsegs[12:132]
        trueSegs[130:140] = flatsegs[133:143]
        self._ui.label_meanSeg.setText("Mean: %.2f" % trueSegs.mean())
        self._ui.label_maxSeg.setText("Maximum: %.2f" % trueSegs.max())
        self._ui.label_minSeg.setText("Minimum: %.2f" % trueSegs.min())



    def shutDown(self):
        pass
        #
        #if self._scanner:
        #    self._scanner.wait()


class ApplyToMirror(QtCore.QThread):
    def __init__(self, control):
        QtCore.QThread.__init__(self)
        self._control = control

    def run(self):
        self._control.applyToMirror()

class ApplyManyZernsToMirror(QtCore.QThread):

    def __init__(self, control, mode, maxAmp, minAmp, num,
                 wTime, radius, mask, clearFirst, withRunning):
        '''
        Applies range of amplitudes of a zernike mode to the mirror

        :Parameters:
            *control*: class
                Deformable mirror control class
            *mode*: int
                Zernike mode to apply
            *maxAmp*: float
                Maximum amplitude of Zernike mode
            *minAmp*: float
                Minimum amplitude of Zernike mode
            *num*: int
                Number of amplitudes to apply to mirror
            *wTime*: float
                Wait time in milliseconds
            *mask*: boolean
                Whether or not to restrict Zernike calculated pattern
                to NA-limited circle or to cover entire square mirror
                       

        '''
        QtCore.QThread.__init__(self)   
        self._control = control
        self.mode = mode
        self.maxAmp = maxAmp
        self.minAmp = minAmp
        self.num = num
        self.wTime = wTime
        self.mask = mask
        self.withRunning = withRunning
        self.clearFirst = clearFirst
        self.radius = radius

    def run(self):
        self._control.varyZernAmp(self.mode, self.maxAmp, self.minAmp,
                                  self.num, self.wTime, useMask=self.mask,
                                  radius = self.radius,
                                  clearfirst = self.clearFirst,
                                  externallyCalled = self.withRunning)

class ApplyManyZernRadiiToMirror(QtCore.QThread):
    def __init__(self, control, mode, amp, maxR, minR, num,
                 wTime, radius, mask, clearFirst, withRunning):

        QtCore.QThread.__init__(self)   
        self._control = control
        self.mode = mode
        self.amp = amp
        self.maxR = maxR
        self.minR = minR
        self.num = num
        self.wTime = wTime
        self.mask = mask
        self.withRunning = withRunning
        self.clearFirst = clearFirst
        self.radius = radius

    def run(self):
        self._control.varyZernRadii(self.mode, self.amp, self.maxR, self.minR,
                                    self.num, self.wTime, useMask=self.mask,
                                    radius = self.radius,
                                    clearfirst = self.clearFirst,
                                    externallyCalled = self.withRunning)
        
class ApplyManyOffsetsToGroup(QtCore.QThread):
    def __init__(self, control, maxAmp, minAmp, num, wTime, clearFirst, withRunning):
        QtCore.QThread.__init__(self)
        self._control = control
        self.maxAmp = maxAmp
        self.minAmp = minAmp
        self.num = num
        self.wTime = wTime
        self.withRunning = withRunning
        self.clearFirst = clearFirst

    def run(self):
        self._control.varyGroupOffset(self.maxAmp,self.minAmp,self.num,self.wTime,
                                      clearfirst = self.clearFirst,
                                      externallyCalled = self.withRunning)

class ApplyManyFromFile(QtCore.QThread):
    def __init__(self, control, filename, wTime, withRunning):
        QtCore.QThread.__init__(self)
        self._control = control
        self.wTime = wTime
        self.withRunning = withRunning
        self.filename = filename

    def run(self):
        self._control.varyArbitrary(self.filename, self.wTime,
                                    externallyCalled = self.withRunning)
        

class ApplyManyMultsToMirror(QtCore.QThread):
    def __init__(self, control, minMult, maxMult, num, wTime, withRunning):
        QtCore.QThread.__init__(self)
        self._control = control
        self.wTime = wTime
        self.withRunning = withRunning
        self.minMult = minMult
        self.maxMult = maxMult
        self.num = num

    def run(self):
        self._control.varyMultiplierCurrent(self.minMult, self.maxMult, self.num,
                                            self.wTime, externallyCalled = self.withRunning)

class Modulation:
    def __init__(self, index, ui):
        self.index = index
        self.checkbox = QtGui.QCheckBox(str(self.index))
        self.checkbox.stateChanged.connect(ui._modulation_toggled)
        self.checkbox.toggle()



class FitResultsDialog(QtGui.QDialog):
    
    def __init__(self, PF, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.PF = PF
        self.ui = fit_results_design.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.lineEditCoefficients.setText(str(PF.zernike_coefficients))
        self.ui.mplwidget.figure.delaxes(self.ui.mplwidget.figure.axes[0])
        axes_raw = self.ui.mplwidget.figure.add_subplot(131)
        axes_raw.matshow(PF.phase, cmap='RdBu')
        axes_raw.set_title('Raw data')
        axes_fit = self.ui.mplwidget.figure.add_subplot(132)
        axes_fit.matshow(PF.zernike, cmap='RdBu', vmin=PF.phase.min(), vmax=PF.phase.max())
        axes_fit.set_title('Fit')
        axes_coefficients = self.ui.mplwidget.figure.add_subplot(133)
        axes_coefficients.bar(np.arange(15), PF.zernike_coefficients)
        axes_coefficients.set_title('Zernike coefficients')


    def getRemove(self):
        return self.ui.checkBoxRemovePTTD.isChecked()


class Scanner(QtCore.QThread):

    def __init__(self, control, range_, nSlices, nFrames, center_xy, fname, maskRadius, maskCenter):
        QtCore.QThread.__init__(self)
        
        self.control = control
        self.range_ = range_
        self.nSlices = nSlices
        self.nFrames = nFrames
        self.center_xy = center_xy
        self.fname = fname
        self.maskRadius = maskRadius
        self.maskCenter = maskCenter

    def run(self):
        self.control.acquirePSF(self.range_, self.nSlices, self.nFrames,
                                self.center_xy, self.fname,
                                self.maskRadius, self.maskCenter)

'''
class RunningSharpness(QtCore.QThread):
                                               
    def __init__(self, control, im_size, maskRadius, maskCenter, pixelSize, diffLimit, plotWidget,
                 varyZern=False):
        QtCore.QThread.__init__(self)
        self.control = control
        self.pixelSize = pixelSize
        self.diffLimit = diffLimit
        self.maskRadius = maskRadius
        self.maskCenter = maskCenter
        self.plotWidget = plotWidget

        self._on = True
        self._varyZern = varyZern
        self.numZernsToVary = 100
        self._zerns = 0
                                               

        nx,ny = im_size

        if maskCenter[0] > -1 and maskCenter[1]>-1:
            x,y = np.meshgrid(np.arange(nx),np.arange(ny))
            x -= maskCenter[0]
            y -= maskCenter[1]
            r_pxl = _msqrt(x**2 + y**2)
            mask = r_pxl<maskRadius
            self.mask = mask
        else:
            self.mask = None

    def turnOff(self):
        self._on = False

    def run(self):
        self._zerns = 0
        while self._on:
            sharpness, sharpnessList = self.control.findSharpnessEachFrame(self.pixelSize, self.diffLimit, self.mask)
            if sharpness is not None:
                if not self._varyZern:
                    self.plotWidget.figure.axes[0].plot(sharpnessList)
                    self.plotWidget.draw()
                if self._varyZern:
                    self.emit(QtCore.SIGNAL('nextModulation(int)'),self._zerns)
                    time.sleep(0.1)
                    self._zerns += 1
                    if self._zerns == 101:
                        self.plotWidget.figure.axes[0].plot(sharpnessList[-100:], '--ro')
                        self.plotWidget.draw()
                        self._on = False
                        self.emit(QtCore.SIGNAL('doneAdvancingZern'))

'''
