#!/usr/bin/python


from PyQt4 import QtCore,QtGui
import inLib


class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        print 'Initializing sliScan UI.'
        inLib.ModuleUI.__init__(self, control, ui_control, 'modules.sliscan.sliscan_design')
        self._ui.doubleSpinBoxStart.valueChanged.connect(self.calcScanParams)
        self._ui.doubleSpinBoxEnd.valueChanged.connect(self.calcScanParams)
        self._ui.spinBoxNSteps.valueChanged.connect(self.calcScanParams)
        self._window.connect(self._ui.pushButtonScan,QtCore.SIGNAL('clicked()'), self.scan)
        self._ui.pushButtonCalculateSLI.clicked.connect(self.calculateScanPatterns)

        self._ui.labelFile.dragEnterEvent = self._drag_enter_event
        self._ui.labelFile.dragMoveEvent = self._drag_move_event
        self._ui.labelFile.dropEvent = self._drop_event

        self.calcScanParams()

        self._calculator = None
        self._scanner = None

        print 'sliScan UI initialized.'


    def _drag_enter_event(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()


    def _drag_move_event(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()


    def _drop_event(self, event):
        for url in event.mimeData().urls():
            path = str(url.toLocalFile())
            self._ui.labelFile.setText(path)


    def getParams(self):
        start = self._ui.doubleSpinBoxStart.value()
        end = self._ui.doubleSpinBoxEnd.value()
        nSteps = self._ui.spinBoxNSteps.value()
        nFrames = self._ui.spinBoxNFrames.value()
        return start, end, nSteps, nFrames


    def calcScanParams(self):
        start, end, nSteps, nFrames = self.getParams()
        step = self._control.calcScanParams(start, end, nSteps)
        self._ui.labelStepSize.setText(str(step))


    def calculateScanPatterns(self):
        if not self._calculator or not self._calculator.isRunning():
            filename = str(QtGui.QFileDialog.getSaveFileName(None,'Save to file','','*.npy'))
            if filename:
                start, end, nSteps, nFrames = self.getParams()
                self._calculator = Calculator(self._control, start, end, nSteps, filename)
                self._calculator.start()


    def scan(self):
        if not self._scanner or not self._scanner.isRunning():
            scan_file = str(self._ui.labelFile.text())
            if scan_file != '-':
                start, end, nSteps, nFrames = self.getParams()
                if self._ui.checkBoxSave.isChecked():
                    filename = str(self._ui.lineEditSave.text())
                else:
                    filename = None
                self._ui.pushButtonScan.setEnabled(False)
                self._scanner = Scanner(self._control, scan_file, nFrames, filename)
                self._window.connect(self._scanner, QtCore.SIGNAL('sliscanDone'),
                        self._on_scan_done)
                self._scanner.start()


    def _on_scan_done(self):
        self._ui.pushButtonScan.setEnabled(True)


    def shutDown(self):
        if self._calculator:
            self._calculator.wait()
        if self._scanner:
            self._scanner.wait()



class Calculator(QtCore.QThread):

    def __init__(self, control, start, end, nSteps, filename):
        QtCore.QThread.__init__(self)

        self.control = control
        self.start_ = start
        self.end = end
        self.nSteps = nSteps
        self.filename = filename


    def run(self):
        self.control.calculateSLIScan(self.start_, self.end, self.nSteps, self.filename)



class Scanner(QtCore.QThread):

    def __init__(self, control, scan_file, nFrames, filename):
        QtCore.QThread.__init__(self)

        self.control = control
        self.scan_file = scan_file
        self.nFrames = nFrames
        self.filename = filename


    def run(self):
        self.control.scan(self.scan_file, self.nFrames, self.filename)
        self.emit(QtCore.SIGNAL('sliscanDone'))
