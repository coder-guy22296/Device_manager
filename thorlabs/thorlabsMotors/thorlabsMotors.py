#!/usr/bin/python
#
import inLib
import sys
import time

class Control(inLib.Device):
    def __init__(self, settings):
        inLib.Device.__init__(self, 'thorlabs.thorlabsMotors.thorlabsMotors_api', settings)

    def initServo(self):
         self._api.initServo()

    def jogUp(self):
        self._api.jogUp()

    def jogDown(self):
        self._api.jogDown()


    def setStage(self, pos):
        self._api.setStage(pos)

    def shutDown(self):
        self._api.shutDown()

