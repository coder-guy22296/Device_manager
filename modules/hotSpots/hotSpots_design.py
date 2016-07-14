# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'modules\hotSpots\hotSpots_design.ui'
#
# Created: Wed May 23 17:52:41 2012
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(180, 256)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "HotSpots", None, QtGui.QApplication.UnicodeUTF8))
        self.gridLayout_2 = QtGui.QGridLayout(Form)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.listWidgetHotSpots = QtGui.QListWidget(Form)
        self.listWidgetHotSpots.setObjectName(_fromUtf8("listWidgetHotSpots"))
        self.gridLayout_2.addWidget(self.listWidgetHotSpots, 0, 0, 1, 3)
        self.pushButtonAdd = QtGui.QPushButton(Form)
        self.pushButtonAdd.setText(QtGui.QApplication.translate("Form", "Add", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonAdd.setObjectName(_fromUtf8("pushButtonAdd"))
        self.gridLayout_2.addWidget(self.pushButtonAdd, 1, 1, 1, 1)
        self.pushButtonRemove = QtGui.QPushButton(Form)
        self.pushButtonRemove.setText(QtGui.QApplication.translate("Form", "Remove", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonRemove.setObjectName(_fromUtf8("pushButtonRemove"))
        self.gridLayout_2.addWidget(self.pushButtonRemove, 1, 2, 1, 1)
        spacerItem = QtGui.QSpacerItem(0, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem, 1, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

