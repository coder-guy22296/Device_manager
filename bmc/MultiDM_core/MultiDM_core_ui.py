#!/usr/bin/python


from PyQt4 import QtGui,QtCore
import inLib
from Utilities import QExtensions as qext
import numpy as np
from numpy.lib.scimath import sqrt as _msqrt
import copy
import time
import skimage

"""
The required buttons:
    getGeometry
    
    
"""


class UI(inLib.DeviceUI):
    
    def __init__(self, control):

        #path to design_path
        design_path = 'bmc.MultiDM_core.mirror_design_core'
        inLib.DeviceUI.__init__(self, control, design_path)

        self._ui.pushButton_load.clicked.connect(self.loadPattern)
        self._ui.pushButton_rot90.clicked.connect(self.patternRot90)
        self._ui.pushButton_fliplr.clicked.connect(self.patternFlipLR)
        self._ui.pushButton_flipud.clicked.connect(self.patternFlipUD)
        self._ui.pushButton_rotate.clicked.connect(self.patternRotate)
        self._ui.pushButton_getSegs.clicked.connect(self.getSegments)
        self._ui.pushButton_toMirror.clicked.connect(self.toMirror)
        self._ui.pushButton_reconfig.clicked.connect(self.reconfig)
        self._ui.pushButton_mult.clicked.connect(self.setMultiplier)
        self._ui.pushButton_premult.clicked.connect(self.setPreMultiplier)
        self._ui.pushButton_poke.clicked.connect(self.pokeSegment)
        self._ui.pushButton_clear.clicked.connect(self.clearPattern)
        self._ui.pushButton_refresh.clicked.connect(self.refreshPattern)
        self._ui.pushButton_pad.clicked.connect(self.padZeros)
        self._ui.pushButton_applyZern.clicked.connect(self.calcZernike)
        
        # This is replaced by add_zernike to modulate
#         self._ui.pushButton_modulateZernike.clicked.connect(self.modZernike)
        self._ui.pushButton_createGroup.clicked.connect(self.createGroup)
        self._ui.pushButton_setToGroup.clicked.connect(self.setGroupVal)


        self._ui.lineEdit_loadMult.setText("10")
        self._ui.lineEdit_npixels.setText(str(self._control.pixels))
        self._ui.lineEdit_zernAmp.setText("0")
        self._ui.lineEdit_premult.setText(str(self._control.preMultiplier))

        self.pattern=None
        self._modulations = [] # added by Dan on 07/14 
        
        self._ui.pushButton_setMods.clicked.connect(self.set_modulations) # added by Dan
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

    # 07/14: How is the pattern passed from adaptive optics to DM?
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

    def set_modulations(self):
#        07/14: This should serve as "set" in the adaptive optics module
#        adapted from adaptiveOptics_ui
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setMod_status(m.index, state)
            
        self._control.mod_from_pool()
        self.refreshPattern()
        
#             
        # 07/17: complete the modulation
   

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

    
    def modulate(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
#        self._control.modulatePF(self.use_zernike)
    

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
