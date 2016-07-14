# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'scripts_design.ui'
#
# Created: Fri May 25 09:31:26 2012
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
        Form.resize(601, 552)
        Form.setAcceptDrops(False)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Scripts", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.listWidgetFiles = QtGui.QListWidget(Form)
        self.listWidgetFiles.setAcceptDrops(True)
        self.listWidgetFiles.setDragEnabled(False)
        self.listWidgetFiles.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.listWidgetFiles.setObjectName(_fromUtf8("listWidgetFiles"))
        self.gridLayout.addWidget(self.listWidgetFiles, 1, 0, 1, 1)
        self.label = QtGui.QLabel(Form)
        self.label.setText(QtGui.QApplication.translate("Form", "Files:", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.pushButtonLoad = QtGui.QPushButton(Form)
        self.pushButtonLoad.setText(QtGui.QApplication.translate("Form", "Load", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonLoad.setObjectName(_fromUtf8("pushButtonLoad"))
        self.horizontalLayout_2.addWidget(self.pushButtonLoad)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.gridLayout.addLayout(self.horizontalLayout_2, 2, 0, 1, 1)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setText(QtGui.QApplication.translate("Form", "Script:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 1, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.pushButtonNew = QtGui.QPushButton(Form)
        self.pushButtonNew.setText(QtGui.QApplication.translate("Form", "New", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonNew.setObjectName(_fromUtf8("pushButtonNew"))
        self.horizontalLayout.addWidget(self.pushButtonNew)
        self.pushButtonSave = QtGui.QPushButton(Form)
        self.pushButtonSave.setText(QtGui.QApplication.translate("Form", "Save", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonSave.setObjectName(_fromUtf8("pushButtonSave"))
        self.horizontalLayout.addWidget(self.pushButtonSave)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.pushButtonRun = QtGui.QPushButton(Form)
        self.pushButtonRun.setText(QtGui.QApplication.translate("Form", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonRun.setObjectName(_fromUtf8("pushButtonRun"))
        self.horizontalLayout.addWidget(self.pushButtonRun)
        self.gridLayout.addLayout(self.horizontalLayout, 2, 1, 1, 1)
        self.gridLayout.setColumnStretch(0, 3)
        self.gridLayout.setColumnStretch(1, 10)
        self.verticalLayout.addLayout(self.gridLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

