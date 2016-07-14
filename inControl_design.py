# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'inControl_design.ui'
#
# Created: Wed Aug 01 11:30:34 2012
#      by: PyQt4 UI code generator 4.8.6
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
        Form.resize(197, 355)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(Form)
        self.label.setText(QtGui.QApplication.translate("Form", "Devices:", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.listWidgetDevices = QtGui.QListWidget(Form)
        self.listWidgetDevices.setObjectName(_fromUtf8("listWidgetDevices"))
        self.verticalLayout.addWidget(self.listWidgetDevices)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setText(QtGui.QApplication.translate("Form", "Modules:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.listWidgetModules = QtGui.QListWidget(Form)
        self.listWidgetModules.setObjectName(_fromUtf8("listWidgetModules"))
        self.verticalLayout.addWidget(self.listWidgetModules)
        self.label_3 = QtGui.QLabel(Form)
        self.label_3.setText(QtGui.QApplication.translate("Form", "Working Directory:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.labelDir = QtGui.QLabel(Form)
        self.labelDir.setText(QtGui.QApplication.translate("Form", "-", None, QtGui.QApplication.UnicodeUTF8))
        self.labelDir.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.labelDir.setObjectName(_fromUtf8("labelDir"))
        self.horizontalLayout.addWidget(self.labelDir)
        self.pushButtonDir = QtGui.QPushButton(Form)
        self.pushButtonDir.setMaximumSize(QtCore.QSize(25, 16777215))
        self.pushButtonDir.setText(QtGui.QApplication.translate("Form", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonDir.setObjectName(_fromUtf8("pushButtonDir"))
        self.horizontalLayout.addWidget(self.pushButtonDir)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.widget = QtGui.QWidget(Form)
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.lineEdit_workingDir = QtGui.QLineEdit(self.widget)
        self.lineEdit_workingDir.setObjectName(_fromUtf8("lineEdit_workingDir"))
        self.horizontalLayout_2.addWidget(self.lineEdit_workingDir)
        self.pushButton_setDir = QtGui.QPushButton(self.widget)
        self.pushButton_setDir.setMaximumSize(QtCore.QSize(30, 50))
        self.pushButton_setDir.setText(QtGui.QApplication.translate("Form", "Set", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_setDir.setObjectName(_fromUtf8("pushButton_setDir"))
        self.horizontalLayout_2.addWidget(self.pushButton_setDir)
        self.verticalLayout.addWidget(self.widget)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

