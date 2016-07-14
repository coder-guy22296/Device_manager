#!/usr/bin/python


from PyQt4 import QtGui,QtCore
import inLib
from Utilities import QExtensions as qext
import numpy as np
from numpy.lib.scimath import sqrt as _msqrt
import copy
import time

class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):

        #path to design_path
        design_path = 'modules.adaptiveOptics.adaptiveOptics_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        self._ui.doubleSpinBoxRange.setValue(self._control._settings['range'])
        self._ui.spinBoxSlices.setValue(self._control._settings['nSlices'])
        self._ui.spinBoxFrames.setValue(self._control._settings['nFrames'])
        if self._control._settings['filename']:
            self._ui.checkBoxSave.setCheckState(2)
            self._ui.lineEditFile.setText(self._control._settings['filename'])

        self._ui.spinBoxIterations.setValue(self._control._settings['nIterations'])

        self._window.connect(self._ui.pushButtonPSF,QtCore.SIGNAL('clicked()'),self.acquirePSF)
        self._window.connect(self._ui.pushButtonModulate,QtCore.SIGNAL('clicked()'),self.modulate)
        self._window.connect(self._ui.pushButtonSave,QtCore.SIGNAL('clicked()'),self.savePF)

        self._ui.pushButton_modUnwrapped.clicked.connect(self.modulateUnwrapped)

        self._ui.groupBoxModulations.toggled.connect(self._modulations_toggled)

        self._ui.pushButton_stopSharpness.clicked.connect(self.stopSharpness)
        self._ui.pushButton_runningSharpness.clicked.connect(self.runningSharpness)
        self._ui.pushButton_sharpnessVsZern.clicked.connect(self.sharpnessVsZern)

        self._ui.pushButton_unwrap.clicked.connect(self.unwrap)

        self._ui.pushButton_zernFitUnwrapped.clicked.connect(self.fitUnwrapped)
        self._ui.pushButton_modulateZernike.clicked.connect(self.modulateUnwrappedZernike)

        self._modulations = []
        self.use_zernike = False
        self.remove_PTTD = True

        self._ui.mplwidgetPhase.figure.axes[0].hold(False)
        self._ui.mplwidgetPhase_2.figure.axes[0].hold(False)

        self._scanner = None

        self._runningSharpness = np.array([])
        self._sharpnessPlot = None

        self.imsize = (256,256)
        self.pixelSize = 163
        self.diffLimit = 800

        zmin,zmax = self._control.getZernMinMax()
        self._ui.doubleSpinBox_zernAmpMin.setValue(zmin)
        self._ui.doubleSpinBox_zernAmpMax.setValue(zmax)

        self._ui.tabWidgetPF.setEnabled(True)
        self._ui.pushButton_modulateZernike.setEnabled(False)

    def _displayPhase(self, phase):
        self._ui.mplwidgetPhase.figure.axes[0].matshow(phase, cmap='RdBu')
        self._ui.mplwidgetPhase.draw()

    def _plotSharpness(self, sharpness):
        self._ui.mplwidgetSharpness.figure.axes[0].plot(sharpness)
        self._ui.mplwidgetSharpness.draw()


    def _modulation_toggled(self, state):
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setModulationActive(m.index, state)
        self._ui_control.slm.updateModulationDisplay()


    def _modulations_toggled(self, state):
        if state == False:
            for m in self._modulations:
                self._control.setModulationActive(m.index, state)
        else:
            self._modulation_toggled(state)
        self._ui_control.slm.updateModulationDisplay()


    def acquirePSF(self):
        range_ = self._ui.doubleSpinBoxRange.value()
        nSlices = self._ui.spinBoxSlices.value()
        nFrames = self._ui.spinBoxFrames.value()
        center_xy = self._ui.checkBoxCenterLateral.isChecked()
        maskRadius = self._ui.spinBox_maskRadius.value()
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        save = self._ui.checkBoxSave.isChecked()
        fname = None
        if save:
            fname = str(self._ui.lineEditFile.text())
        self._scanner = Scanner(self._control, range_, nSlices, nFrames, center_xy,
                                fname, maskRadius, (cX,cY))
        self._scanner.finished.connect(self._on_scan_done)
        self._ui.pushButtonPSF.setEnabled(False)
        self._scanner.start()


    def _on_scan_done(self):
        self._ui.pushButtonPSF.setEnabled(True)
        self._ui.groupBoxPhase.setEnabled(True)
        time.sleep(2)
        sharpness = self._control.getSharpness()
        self._plotSharpness(sharpness)

    def runningSharpness(self):
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        maskRadius = self._ui.spinBox_maskRadius.value()
        self._sharpnessPlot = RunningSharpness(self._control, self.imsize, maskRadius,
                                               (cX,cY), self.pixelSize, self.diffLimit,
                                               self._ui.mplwidgetSharpness2)
        #self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('nextModulation'),self.advanceModulation)
        self._sharpnessPlot.start()

    def sharpnessVsZern(self):
        zmin = self._ui.doubleSpinBox_zernAmpMin.value()
        zmax = self._ui.doubleSpinBox_zernAmpMax.value()
        self._control.setZernAmps(zmin,zmax)
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        maskRadius = self._ui.spinBox_maskRadius.value()
        self._control.resetZernAmpIndex()
        self._sharpnessPlot = RunningSharpness(self._control, self.imsize, maskRadius,
                                               (cX,cY), self.pixelSize, self.diffLimit,
                                               self._ui.mplwidgetSharpness2, varyZern=True)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('nextModulation(int)'),self.advanceModulation)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('doneAdvancingZern'),self.stopSharpness)
        self._sharpnessPlot.start()

    def stopSharpness(self):
        if self._sharpnessPlot is not None:
            self._sharpnessPlot.turnOff()
            self._sharpnessPlot.wait()
            self._sharpnessPlot.quit()
            self._sharpnessPlot = None

    def advanceModulation(self, mod_index):
        coeff = self._control.advanceModulation()
        self._ui_control.slm.updateModulationDisplay()
        self._ui.label_mod_index.setText("Index: %i" % mod_index)
        self._ui.label_mod_value.setText("Value: %.2f" % coeff)

    def fitPF(self):
        PF = self._control.fit()
        fit_result_dialog = FitResultsDialog(PF)
        if fit_result_dialog.exec_():
            print 'adaptiveOptics: Fit accepted.'
            self.use_zernike = True
            remove = fit_result_dialog.getRemove()
            if remove:
                print 'adaptiveOptics: Remove pistion, tip, tilt and defocus.'
                PF = self._control.removePTTD()
            self._displayPhase(PF.zernike)
        else:
            self.use_zernike = False


    def retrievePF(self):
        pxlSize = self._ui.doubleSpinBoxPixel.value()
        l = self._ui.doubleSpinBoxWavelength.value()
        n = self._ui.doubleSpinBoxIndex.value()
        NA = self._ui.doubleSpinBoxNA.value()
        f = self._ui.doubleSpinBoxFocal.value()
        numWaves = self._ui.spinBox_numWavelengths.value()
        neglect_defocus = self._ui.checkBoxNeglectDefocus.isChecked()
        nIt = self._ui.spinBoxIterations.value()
        invertPF = self._ui.checkBox_invertPF.isChecked()
        resetAmp = self._ui.checkBox_resetAmp.isChecked()

        checked = self._ui.buttonGroupGuess.checkedButton()
        if checked == self._ui.radioButtonPlane:
            guess = ('plane',)
        elif checked == self._ui.radioButtonMirror:
            z0 = self._ui.doubleSpinBoxMirrorDistance.value()
            guess = ('mirror',z0)
        elif checked == self._ui.radioButtonFromFile:
            filename = QtGui.QFileDialog.getOpenFileName(None,'Open initial guess',
                                                  '','*.npy')
            if filename:
                guess = ('file', str(filename))
        if (guess[0] != 'file') or (guess[0] == 'file' and filename):
            PF = self._control.retrievePF(pxlSize, l, n, NA, f, guess, nIt, neglect_defocus, invert=invertPF,
                                          wavelengths=numWaves, resetAmp=resetAmp)
            self._ui.tabWidgetPF.setEnabled(True)
            self._ui.pushButton_modulateZernike.setEnabled(False)
            self._displayPhase(PF)
        self.use_zernike = False 


    def modulate(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modulatePF(self.use_zernike)
        self._ui_control.slm.updateModulationDisplay()

    def modulateUnwrapped(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modulatePF_unwrapped()
        self._ui_control.slm.updateModulationDisplay()

    def savePF(self):
        filename = QtGui.QFileDialog.getSaveFileName(None,'Save to file',
                                                  '','*.npy')
        if filename:
            self._control.savePF(str(filename))

    def unwrap(self):
        unwrappedPhase = self._control.unwrap()
        self._displayPhase(unwrappedPhase)

    def fitUnwrapped(self):
        resultFit = self._control.zernFitUnwrapped()
        self._ui.mplwidgetPhase_2.figure.axes[0].matshow(resultFit, cmap='RdBu')
        self._ui.mplwidgetPhase_2.draw()
        self._ui.pushButton_modulateZernike.setEnabled(True)

    def modulateUnwrappedZernike(self):
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modZernFitUnwrapped()
        self._ui_control.slm.updateModulationDisplay()


    def shutDown(self):
        if self._scanner:
            self._scanner.wait()


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

