#!/usr/bin/python

from PyQt4 import QtCore,QtGui,Qwt5
from functools import partial
import inLib
import sys
import numpy as np

from Utilities import QExtensions as qext


class UI(inLib.ModuleUI):
    
    def __init__(self, control, ui_control):
        design_path = 'modules.testing.testing2_design'
        inLib.ModuleUI.__init__(self, control, ui_control, design_path)

        self.cam_window = self._ui_control.camera._ui

        self._ui.pushButton_placeTab.clicked.connect(self.place)
        self._ui.pushButton_delTab.clicked.connect(self.delete)
        self._ui.pushButton_close.clicked.connect(self.close)

        #self.place()
        #self.close()

        self.index=0

    def place(self):
        self.index = self.cam_window.tabWidget.addTab(self._ui.grpBox, "Testing")

    def delete(self):
        if not self.index:
            self.cam_window.tabWidget.removeTab(self.index)

    def close(self):
        self._window.close()
