#!/usr/bin/python


from PyQt4 import QtGui,QtCore
import inLib
from Utilities import QExtensions as qext
import numpy as np
from numpy.lib.scimath import sqrt as _msqrt
import fit_results_design
import copy
import time
from libs import scipy_gaussfitter

class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        design_path = 'modules.adaptiveOptics.adaptiveOptics_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        self._ui.buttonGroupGuess = QtGui.QButtonGroup()
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonPlane)
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonMirror)
        self._ui.buttonGroupGuess.addButton(self._ui.radioButtonFromFile)
        self._ui.pushButtonFit.clicked.connect(self.fitPF)
        self._ui.pushButtonPhase.clicked.connect(self.retrievePF)

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

        self._ui.groupBoxModulations.toggled.connect(self._modulations_toggled)
        self._ui.pushButton_stopSharpness.clicked.connect(self.stopSharpness)
        self._ui.pushButton_runningSharpness.clicked.connect(self.runningSharpness)
        self._ui.pushButton_sharpnessVsZern.clicked.connect(self.sharpnessVsZern)

        self._ui.pushButton_unwrap.clicked.connect(self.unwrap)

        self._ui.pushButton_zernFitUnwrapped.clicked.connect(self.fitUnwrapped)
        self._ui.spinBox_zernModesToFit.setValue(self._control.zernModesToFit)
        self._ui.spinBox_zernModesToFit.valueChanged.connect(self.setZernModesToFit)
        self._ui.lineEdit_diffLimit.returnPressed.connect(self.setDiffLimit)
        
        # 07/14: this is how the modulation is set to the mirror.        
        self._ui.pushButton_setMods.clicked.connect(self.set_modulations)

        self._ui.pushButton_setZernRadius.clicked.connect(self._setZernikeRadius)
        self.zernRadius = 0

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
        self._ui.lineEdit_diffLimit.setText("%i" % self.diffLimit)

        zmin,zmax = self._control.getZernMinMax()
        self._ui.doubleSpinBox_zernAmpMin.setValue(zmin)
        self._ui.doubleSpinBox_zernAmpMax.setValue(zmax)

        self._ui.tabWidgetPF.setEnabled(True)
        self._ui.pushButton_modulateZernike.setEnabled(False)

        self.hasSLM = self._control.hasSLM
        self.hasMirror = self._control.hasMirror

    def _displayPhase(self, phase):
        self._ui.mplwidgetPhase.figure.axes[0].matshow(phase, cmap='RdBu')
        self._ui.mplwidgetPhase.draw()

    def _plotSharpness(self, sharpness):
        self._ui.mplwidgetSharpness.figure.axes[0].plot(sharpness)
        self._ui.mplwidgetSharpness.draw()

    def _updateImSize(self):
        self.imsize = self._control.updateImSize()


    def _modulation_toggled(self, state):
        pass
        '''
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setModulationActive(m.index, state)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()
        '''


    def _modulations_toggled(self, state):
        pass
        '''
        if state == False:
            for m in self._modulations:
                self._control.setModulationActive(m.index, state)
        else:
            self._modulation_toggled(state)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()
        '''

    def set_modulations(self):
        for m in self._modulations:
            state = m.checkbox.isChecked()
            self._control.setModulationActive(m.index, state)

    def _setZernikeRadius(self):
        radius = int(self._ui.lineEdit_zernRadius.text())
        self.zernRadius = radius

    def setDiffLimit(self):
        self.diffLimit = int(self._ui.lineEdit_diffLimit.text())


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
        '''
        Continously run shaprness
        '''
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        maskRadius = self._ui.spinBox_maskRadius.value()
        self._sharpnessPlot = RunningSharpness(self._control, self.imsize, maskRadius,
                                               (cX,cY), self.pixelSize, self.diffLimit,
                                               self._ui.mplwidgetSharpness2)
        #self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('nextModulation'),self.advanceModulation)
        self._sharpnessPlot.start()

    def sharpnessVsZern(self):
        '''
        started on pushButton_sharpnessVsZern.clicked

        Starts thread RunningSharpness.
        Since varyZern=True when calling thread, won't keep calculating sharpness until stopped
        '''
        self._updateImSize()
        zmin = self._ui.doubleSpinBox_zernAmpMin.value() #minimum zernike amplitude
        zmax = self._ui.doubleSpinBox_zernAmpMax.value() #maximum zernike amplitude
        self._control.setZernAmps(zmin,zmax,num=100)
        cX = int(self._ui.lineEdit_cX.text())
        cY = int(self._ui.lineEdit_cY.text())
        wTime = float(self._ui.lineEdit_waitTime.text())
        maskRadius = self._ui.spinBox_maskRadius.value()
        self._control.resetZernAmpIndex()
        #if self.hasMirror:
        #    self._control.mirror.varyZernAmp()
        numToVary = self._control.getNumberOfZernToVary()
        self._sharpnessPlot = RunningSharpness(self._control, self.imsize, maskRadius,
                                               (cX,cY), self.pixelSize, self.diffLimit,
                                               self._ui.mplwidgetSharpness2, numToVary = numToVary, varyZern=True,
                                               wait_time = wTime)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('nextModulation(int)'),self.advanceModulation)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('doneAdvancingZern'),self.stopSharpness)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('maxArgSharpness(int)'),self.foundMaxArgSharp)
        self._window.connect(self._sharpnessPlot,QtCore.SIGNAL('maxFitSharpness(float)'),self.foundMaxFitSharp)
        self._sharpnessPlot.start()

    def stopSharpness(self):
        '''
        Connected to signal 'doneAdvancingZern' which is called by RunningSharpness thread

        Just stops the thread and sets it back to None
        '''
        if self._sharpnessPlot is not None:
            self._sharpnessPlot.turnOff()
            self._sharpnessPlot.wait()
            self._sharpnessPlot.quit()
            self._sharpnessPlot = None

    def foundMaxArgSharp(self, argmax):
        '''
        Connects to signal 'maxArgSharpness(int)'
        '''
        self._ui.label_sharpnessArgMax.setText("Arg. max: %i" % argmax)

    def foundMaxFitSharp(self, fit):
        self._ui.label_sharpnessFitMax.setText("Fit max: %.2f" % fit)

    def advanceModulation(self, mod_index):
        '''
        Connected to signal 'nextModulation(int)' which is called by RunningSharpness

        This is triggered each time a new pattern is to be displayed on the adaptive optics device.
        '''
        coeff = self._control.advanceModulation()
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay() #updates display
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
        symmeterize = self._ui.checkBox_symmeterize.isChecked()

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
                                          wavelengths=numWaves, resetAmp=resetAmp,
                                          symmeterize=symmeterize)
            
#             def retrievePF(self, dx, l, n, NA, f, guess, nIt, neglect_defocus=True,
#                    invert=False, wavelengths=1, resetAmp=True,
#                    symmeterize=False):
            self._ui.tabWidgetPF.setEnabled(True)
            self._ui.pushButton_modulateZernike.setEnabled(False)
            self._displayPhase(PF)
        self.use_zernike = False 

    # This is connection from the button "Modulate"
    # 07/18: Adapt this to the MultiDM_core
    def modulate(self):
        # 07/21: remove modulateUnwrapped, because the phase is by default unwrapped.
        modulation = Modulation(len(self._modulations), self)
        self._ui.verticalLayoutModulations.insertWidget(0, modulation.checkbox)
        self._modulations.append(modulation)
        self._control.modulatePF(self.use_zernike)
        if self.hasSLM:
            self._ui_control.slm.updateModulationDisplay()
            
            

    def savePF(self):
        filename = QtGui.QFileDialog.getSaveFileName(None,'Save to file',
                                                  '','*.npy')
        if filename:
            self._control.savePF(str(filename))

    def unwrap(self):
        unwrappedPhase = self._control.unwrap()
        self._displayPhase(unwrappedPhase)

    def setZernModesToFit(self):
        nmodes = self._ui.spinBox_zernModesToFit.value()
        self._control.setZernModesToFit(nmodes)

    def fitUnwrapped(self):
        ignore4 = self._ui.checkBox_ignore4.isChecked()
        resultFit = self._control.zernFitUnwrapped(skip4orders=ignore4)
        self._ui.mplwidgetPhase_2.figure.axes[0].matshow(resultFit, cmap='RdBu')
        self._ui.mplwidgetPhase_2.draw()
        self._ui.pushButton_modulateZernike.setEnabled(True)


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
                 numToVary = 100, varyZern=False, wait_time = 0.1):
        '''
        Finds sharpness

        Parameters
            :control: class
            :im_size: (int,int)
            :maskRadius: int
            :maskCenter: int
            :pixelSize: int
            :diffLimit: int
                Diffraction limited resolution (in same units as pixelSize)
            :plotWidget:
            :numToVary: int
                Optional
            :varyZern: boolean
                Optional
            :wait_time: float
                Optional. Time between patterns applied to adaptive optics device (seconds)

        Emits signal 'nextModulation(int)' when ready to add new modulation...
            nextModulation takes parameter that counts upward from 0
        Emits signal 'doneAdvancingZern' when done
        

        '''
        QtCore.QThread.__init__(self)
        self.control = control
        self.pixelSize = pixelSize
        self.diffLimit = diffLimit
        self.maskRadius = maskRadius
        self.maskCenter = maskCenter
        self.plotWidget = plotWidget
        self.wait_time = wait_time
        self._on = True
        self._varyZern = varyZern
        self.numZernsToVary = numToVary
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
                    time.sleep(self.wait_time)
                    self._zerns += 1
                    if self._zerns == self.numZernsToVary+1:
                        sharpToPlot = sharpnessList[-1*self.numZernsToVary:]
                        maxSharpness = np.argmax(sharpToPlot)
                        start = maxSharpness-20
                        if start<0:
                            start=0
                        stop = maxSharpness+20
                        if stop>self.numZernsToVary:
                            stop = self.numZernsToVary
                        portion = sharpToPlot[start:stop]
                        gaussfit = scipy_gaussfitter.fitgaussian1d(np.array(portion))
                        g = scipy_gaussfitter.gaussian1d(gaussfit[0],gaussfit[1],gaussfit[2])
                        xaxis = np.arange(start,stop)
                        fitResult = g(xaxis)
                        self.emit(QtCore.SIGNAL('maxFitSharpness(float)'),gaussfit[1]+start)
                        self.emit(QtCore.SIGNAL('maxArgSharpness(int)'),maxSharpness)
                        self.plotWidget.figure.axes[0].plot(sharpToPlot, '--ro')
                        #self.plotWidget.figure.axes[0].plot(xaxis,fitResult,'-k')
                        self.plotWidget.draw()
                        self._on = False
                        self.emit(QtCore.SIGNAL('doneAdvancingZern'))

