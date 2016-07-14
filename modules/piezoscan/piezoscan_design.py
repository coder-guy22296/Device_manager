# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\piezoscan_design.ui'
#
# Created: Wed Apr 18 14:50:08 2012
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
        Form.resize(311, 119)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Piezoscan", None, QtGui.QApplication.UnicodeUTF8))
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_7 = QtGui.QLabel(Form)
        self.label_7.setText(QtGui.QApplication.translate("Form", "Start [um]:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout.addWidget(self.label_7, 3, 0, 1, 1)
        self.label_8 = QtGui.QLabel(Form)
        self.label_8.setText(QtGui.QApplication.translate("Form", "End [um]:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout.addWidget(self.label_8, 4, 0, 1, 1)
        self.label_9 = QtGui.QLabel(Form)
        self.label_9.setText(QtGui.QApplication.translate("Form", "# Steps:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout.addWidget(self.label_9, 5, 0, 1, 1)
        self.doubleSpinBoxStart = QtGui.QDoubleSpinBox(Form)
        self.doubleSpinBoxStart.setMinimum(-99.99)
        self.doubleSpinBoxStart.setProperty("value", -2.0)
        self.doubleSpinBoxStart.setObjectName(_fromUtf8("doubleSpinBoxStart"))
        self.gridLayout.addWidget(self.doubleSpinBoxStart, 3, 1, 1, 1)
        self.doubleSpinBoxEnd = QtGui.QDoubleSpinBox(Form)
        self.doubleSpinBoxEnd.setMinimum(-99.0)
        self.doubleSpinBoxEnd.setProperty("value", 2.0)
        self.doubleSpinBoxEnd.setObjectName(_fromUtf8("doubleSpinBoxEnd"))
        self.gridLayout.addWidget(self.doubleSpinBoxEnd, 4, 1, 1, 1)
        self.spinBoxNSteps = QtGui.QSpinBox(Form)
        self.spinBoxNSteps.setMinimum(2)
        self.spinBoxNSteps.setMaximum(9999)
        self.spinBoxNSteps.setProperty("value", 5)
        self.spinBoxNSteps.setObjectName(_fromUtf8("spinBoxNSteps"))
        self.gridLayout.addWidget(self.spinBoxNSteps, 5, 1, 1, 1)
        self.pushButtonScan = QtGui.QPushButton(Form)
        self.pushButtonScan.setEnabled(True)
        self.pushButtonScan.setText(QtGui.QApplication.translate("Form", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonScan.setObjectName(_fromUtf8("pushButtonScan"))
        self.gridLayout.addWidget(self.pushButtonScan, 6, 3, 1, 1)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setText(QtGui.QApplication.translate("Form", "# Frames/step:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 6, 0, 1, 1)
        self.spinBoxNFrames = QtGui.QSpinBox(Form)
        self.spinBoxNFrames.setMinimum(1)
        self.spinBoxNFrames.setMaximum(100)
        self.spinBoxNFrames.setProperty("value", 1)
        self.spinBoxNFrames.setObjectName(_fromUtf8("spinBoxNFrames"))
        self.gridLayout.addWidget(self.spinBoxNFrames, 6, 1, 1, 1)
        self.label_11 = QtGui.QLabel(Form)
        self.label_11.setText(QtGui.QApplication.translate("Form", "Step size [nm]:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout.addWidget(self.label_11, 3, 2, 1, 1)
        self.labelStepSize = QtGui.QLabel(Form)
        self.labelStepSize.setText(QtGui.QApplication.translate("Form", "-", None, QtGui.QApplication.UnicodeUTF8))
        self.labelStepSize.setObjectName(_fromUtf8("labelStepSize"))
        self.gridLayout.addWidget(self.labelStepSize, 3, 3, 1, 1)
        self.checkBoxSave = QtGui.QCheckBox(Form)
        self.checkBoxSave.setText(QtGui.QApplication.translate("Form", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBoxSave.setChecked(True)
        self.checkBoxSave.setObjectName(_fromUtf8("checkBoxSave"))
        self.gridLayout.addWidget(self.checkBoxSave, 4, 2, 1, 1)
        self.lineEditFile = QtGui.QLineEdit(Form)
        self.lineEditFile.setObjectName(_fromUtf8("lineEditFile"))
        self.gridLayout.addWidget(self.lineEditFile, 4, 3, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

