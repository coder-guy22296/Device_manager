# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'testing2_design.ui'
#
# Created: Fri Aug 31 17:33:17 2012
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
        Form.resize(224, 225)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.centralwidget = QtGui.QWidget(Form)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.pushButton_placeTab = QtGui.QPushButton(self.centralwidget)
        self.pushButton_placeTab.setGeometry(QtCore.QRect(10, 20, 75, 23))
        self.pushButton_placeTab.setText(QtGui.QApplication.translate("Form", "place", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_placeTab.setObjectName(_fromUtf8("pushButton_placeTab"))
        self.pushButton_delTab = QtGui.QPushButton(self.centralwidget)
        self.pushButton_delTab.setGeometry(QtCore.QRect(10, 50, 75, 23))
        self.pushButton_delTab.setText(QtGui.QApplication.translate("Form", "del", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_delTab.setObjectName(_fromUtf8("pushButton_delTab"))
        self.pushButton_close = QtGui.QPushButton(self.centralwidget)
        self.pushButton_close.setGeometry(QtCore.QRect(10, 80, 75, 23))
        self.pushButton_close.setText(QtGui.QApplication.translate("Form", "close", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_close.setObjectName(_fromUtf8("pushButton_close"))
        Form.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(Form)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 224, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setTitle(QtGui.QApplication.translate("Form", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        Form.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(Form)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        Form.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(Form)
        self.toolBar.setWindowTitle(QtGui.QApplication.translate("Form", "toolBar", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        Form.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

