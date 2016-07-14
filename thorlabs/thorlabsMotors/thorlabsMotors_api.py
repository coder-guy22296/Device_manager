#!/usr/bin/python
#
# Active X based APTUser control.
#


from PyQt4 import QtCore, QtGui, QAxContainer

import sys
import time
import mgmotorAX

class API():
    def __init__(self):
        self.stage = mgmotorAX.APTUser1()
        #self.servo = mgmotorAX.APTUserTest
        #self.stage.show()
        self.servo = None

    def initServo(self):
        self.servo = mgmotorAX.APTUserTest()
        self.servo.show()

    def setStage(self, pos):
        pass
        #self.stage.moveTo(pos)

    def jogUp(self):
        self.servo.jogUp()

    def jogDown(self):
        self.servo.jogDown()

    def shutDown(self):
        pass
        #self.stage.close()

