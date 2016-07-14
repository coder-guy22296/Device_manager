#!/usr/bin/python


from PyQt4 import QtCore,QtGui
import inLib



class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        print 'Initializing Piezoscan UI.'
        inLib.ModuleUI.__init__(self, control, ui_control, 'modules.piezoscan.piezoscan_design')
        print 'Piezoscan design initialized.'
        self._ui.doubleSpinBoxStart.valueChanged.connect(self.calcScanParams)
        self._ui.doubleSpinBoxEnd.valueChanged.connect(self.calcScanParams)
        self._ui.spinBoxNSteps.valueChanged.connect(self.calcScanParams)
        self._ui.spinBoxNFrames.valueChanged.connect(self.calcScanParams)
        self._window.connect(self._ui.pushButtonScan,QtCore.SIGNAL('clicked()'),self.scan)
        print 'Calculating piezoscan parameters.'
        self.calcScanParams()
        
        self._scanner = None

        print 'Piezoscan UI initialized.'


    def getParams(self):
        start = self._ui.doubleSpinBoxStart.value()
        end = self._ui.doubleSpinBoxEnd.value()
        nSteps = self._ui.spinBoxNSteps.value()
        nFrames = self._ui.spinBoxNFrames.value()
        return start,end,nSteps,nFrames


    def calcScanParams(self):
        start,end,nSteps,nFrames = self.getParams()
        step = self._control.calcScanParams(start,end,nSteps)
        self._ui.labelStepSize.setText(str(step))
        self._ui.lineEditFile.setText('piezoscan_{0}_{1}_{2}_{3}'.format(start,
            end,nSteps,nFrames))


    def scan(self):
        if self._ui.pushButtonScan.text() == 'Start':
            start,end,nSteps,nFrames = self.getParams()
            save = self._ui.checkBoxSave.isChecked()
            if save:
                filename = str(self._ui.lineEditFile.text())
            else:
                filename = None
            self._scanner = Scanner(self._control, start, end, nSteps, nFrames, filename)
            self._scanner.finished.connect(self._on_done)
            self._ui.pushButtonScan.setText('Stop')
            self._scanner.start()
        else:
            #self._scanner.stop()
            self._control.active = False
            self._scanner.wait()
            self._ui.pushButtonScan.setText('Start')


    def _on_done(self):
        self._ui.pushButtonScan.setText('Start')


    def shutDown(self):
        if self._scanner:
            self._scanner.wait()



class Scanner(QtCore.QThread):

    def __init__(self, control, start, end, nSteps, nFrames, filename):
        QtCore.QThread.__init__(self)
        
        self.control = control
        self.start_ = start
        self.end = end
        self.nSteps = nSteps
        self.nFrames = nFrames
        self.filename = filename


    def run(self):
        self.control.scan(self.start_, self.end, self.nSteps, self.nFrames,
                                           self.filename)

    def stop(self):
        print 'trying to stop'
        self.control.stop()
