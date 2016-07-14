#!/usr/bin/python

from PyQt4 import QtCore, QtGui, QAxContainer
import inLib
import time


class UI(inLib.DeviceUI):
    
    def __init__(self, control):
        design_path = 'thorlabs.thorlabsMotors.thorlabsMotors_design'
        print 'Thorlabs stage: Initializing UI.'
        inLib.DeviceUI.__init__(self, control, design_path)

        #self.homeX, self.homeY, z = self._control.position()
        self.sup = 61
        self.sdown = 31
        self._ui.lineEdit_pos.returnPressed.connect(self._move)
        self._ui.lineEdit_scaninit.returnPressed.connect(self._scaninit)
        self._ui.pushButtonUp.clicked.connect(self._stepUp)
        self._ui.pushButtonDown.clicked.connect(self._stepDown)
        self._ui.pushButtonBL.clicked.connect(self._BLcorrect)
        self.scan_0 = 5.00
        #self._ui.pushButtonHomeXY.clicked.connect(self._moveHomeXY)

        #self._ui.updateInfo_checkBox.stateChanged.connect(self._updateInfo)


        #self._positionUpdater = QtCore.QTimer()
        #self._positionUpdater.timeout.connect(self.updatePositionLabels)
        #self._positionUpdater.start(500) #only if updateInfo checkbox is true

        #self._ui.go_pushButton.clicked.connect(self._sawtooth)

        #self.aptUser = APTUser1()

        #self.aptUser = self._control.servo()

        #self.aptUser.show()

        self._control.initServo()

    def _updateInfo(self):
        state = self._ui.updateInfo_checkBox.isChecked()
        if state:
            self._positionUpdater.start(500)
        else:
            if self._positionUpdater.isActive():
                self._positionUpdater.stop()

    '''
    def _sawtooth(self):
        if self._ui.go_pushButton.text() == 'Go':
            movement = float(self._ui.sawtoothMovement_lineEdit.text())
            repeats = int(self._ui.sawtoothRepeats_lineEdit.text())
            timedelay = int(self._ui.sawtoothTimeDelay_lineEdit.text())
            if self._ui.xsawtooth_radioButton.isChecked():
                xy = 'x'
            else:
                xy = 'y'
            self._sawtoothThread = Sawtooth(self._control, xy,
                                            movement, repeats, timedelay)
            self._window.connect(self._sawtoothThread,
                                 QtCore.SIGNAL('sawtoothDone'),
                                 self._on_sawtoothDone)
            self._ui.go_pushButton.setText('Wait...')
            self._sawtoothThread.start()
    '''
    
        
    def _on_sawtoothDone(self):
        self._ui.go_pushButton.setText('Go')
        self._sawtoothThread = None

    def _move(self):
        pos = float(self._ui.lineEdit_pos.text())
        print pos
        self._control.setStage(pos)
#        self.aptUser.moveTo(pos)
    def _scaninit(self):
        self.scan_0 = float(self._ui.lineEdit_scaninit.text())
        print("Scan starts at: %4.2f " %self.scan_0)

    def _stepUp(self):
        #self.aptUser.jogUp()
        self._control.jogUp()

    def _stepDown(self):
        #self.aptUser.jogDown()
        self._control.jogDown()
        
    def _BLcorrect(self):
#       
        self._control.setStage(self.scan_0+0.0090)
        time.sleep(3)
        for ii in xrange(self.sdown+1):
            self._control.jogDown()
            time.sleep(0.27)

    '''
    def _moveHomeXY(self):
        self._control.goAbsolute(self.homeX,self.homeY)


    def updatePositionLabels(self):
        x,y,z = self._control.position()
        self._ui.labelX.setText('%.3f' % x)
        self._ui.labelY.setText('%.3f' % y)
    '''

    def shutDown(self):
        pass
        #self._positionUpdater.stop()

class APTUser1(QAxContainer.QAxWidget):
    def __init__(self, parent = None):
        QAxContainer.QAxWidget.__init__(self, parent)
        self.setControl("MGMOTOR.MGMotorCtrl.1")
        self.dynamicCall('SetHWSerialNum(int)', 83854883)
        self.dynamicCall('StartCtrl()')

        try:
            self.dynamicCall('LoadParamSet(char[])', "smallstep")
        except:
            print "loadparamset for APT User failed"

    def moveTo(self, pos):
        self.dynamicCall('MoveAbsoluteEx(int, double, double, bool)', 0, pos, 0.0, False)

    def jogUp(self):
        self.dynamicCall('MoveJog(int,int)', 0,1)

    def jogDown(self):
        self.dynamicCall('MoveJog(int,int)', 0,2)


'''
class Sawtooth(QtCore.QThread):
    def __init__(self, control, xy, movement, repeats, timedelay):
        QtCore.QThread.__init__(self)

        self.control = control
        self.movement = movement
        self.repeats = repeats
        self.timedelay = timedelay
        self.xy = xy

    def run(self):
        self.control.sawtooth(self.xy,self.movement,self.repeats,self.timedelay)
        self.emit(QtCore.SIGNAL('sawtoothDone'))
'''
