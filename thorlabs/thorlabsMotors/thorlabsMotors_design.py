# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'thorlabsMotors_design.ui'
#
# Created: Tue Mar 17 16:54:54 2015
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
        Form.resize(450, 384)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Thorlabs stage", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox = QtGui.QGroupBox(Form)
        self.groupBox.setTitle(QtGui.QApplication.translate("Form", "Position Control", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setGeometry(QtCore.QRect(18, 69, 26, 16))
        self.label_4.setText(QtGui.QApplication.translate("Form", "Step:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        
        self.pushButtonBL = QtGui.QPushButton(self.groupBox)
        self.pushButtonBL.setGeometry(QtCore.QRect(350, 60, 85, 23))
        self.pushButtonBL.setText(QtGui.QApplication.translate("Form", "BL_correction", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonBL.setObjectName(_fromUtf8("pushButtonBL"))       
        
        self.pushButtonUp = QtGui.QPushButton(self.groupBox)
        self.pushButtonUp.setGeometry(QtCore.QRect(240, 30, 75, 23))
        self.pushButtonUp.setText(QtGui.QApplication.translate("Form", "^", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonUp.setObjectName(_fromUtf8("pushButtonUp"))
        self.pushButtonHomeXY = QtGui.QPushButton(self.groupBox)
        self.pushButtonHomeXY.setGeometry(QtCore.QRect(240, 60, 75, 23))
        self.pushButtonHomeXY.setText(QtGui.QApplication.translate("Form", "Home", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonHomeXY.setObjectName(_fromUtf8("pushButtonHomeXY"))
        self.pushButtonDown = QtGui.QPushButton(self.groupBox)
        self.pushButtonDown.setGeometry(QtCore.QRect(240, 90, 75, 23))
        self.pushButtonDown.setText(QtGui.QApplication.translate("Form", "v", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButtonDown.setObjectName(_fromUtf8("pushButtonDown"))
        self.doubleSpinBoxStep = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBoxStep.setGeometry(QtCore.QRect(50, 70, 49, 20))
        self.doubleSpinBoxStep.setSingleStep(0.1)
        self.doubleSpinBoxStep.setProperty("value", 1.0)
        self.doubleSpinBoxStep.setObjectName(_fromUtf8("doubleSpinBoxStep"))
        self.lineEdit_pos = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_pos.setGeometry(QtCore.QRect(35, 23, 80, 20))
        self.lineEdit_pos.setObjectName(_fromUtf8("lineEdit_pos"))
        self.lineEdit_scaninit = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_scaninit.setGeometry(QtCore.QRect(125, 23, 80, 20))
        self.lineEdit_scaninit.setObjectName(_fromUtf8("lineEdit_scaninit"))
        
        
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_3 = QtGui.QGroupBox(Form)
        self.groupBox_3.setTitle(QtGui.QApplication.translate("Form", "Information", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_3)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_18 = QtGui.QLabel(self.groupBox_3)
        self.label_18.setText(QtGui.QApplication.translate("Form", "Position [um]", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.gridLayout_3.addWidget(self.label_18, 1, 1, 1, 1)
        self.updateInfo_checkBox = QtGui.QCheckBox(self.groupBox_3)
        self.updateInfo_checkBox.setText(QtGui.QApplication.translate("Form", "update info?", None, QtGui.QApplication.UnicodeUTF8))
        self.updateInfo_checkBox.setObjectName(_fromUtf8("updateInfo_checkBox"))
        self.gridLayout_3.addWidget(self.updateInfo_checkBox, 1, 0, 1, 1)
        self.labelPos = QtGui.QLabel(self.groupBox_3)
        self.labelPos.setText(QtGui.QApplication.translate("Form", "-", None, QtGui.QApplication.UnicodeUTF8))
        self.labelPos.setObjectName(_fromUtf8("labelPos"))
        self.gridLayout_3.addWidget(self.labelPos, 1, 2, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_3)
        self.groupBox_2 = QtGui.QGroupBox(Form)
        self.groupBox_2.setTitle(QtGui.QApplication.translate("Form", "Saw-tooth movement", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.xsawtooth_radioButton = QtGui.QRadioButton(self.groupBox_2)
        self.xsawtooth_radioButton.setText(QtGui.QApplication.translate("Form", "x", None, QtGui.QApplication.UnicodeUTF8))
        self.xsawtooth_radioButton.setObjectName(_fromUtf8("xsawtooth_radioButton"))
        self.gridLayout_2.addWidget(self.xsawtooth_radioButton, 0, 1, 1, 1)
        self.ysawtooth_radioButton = QtGui.QRadioButton(self.groupBox_2)
        self.ysawtooth_radioButton.setText(QtGui.QApplication.translate("Form", "y", None, QtGui.QApplication.UnicodeUTF8))
        self.ysawtooth_radioButton.setObjectName(_fromUtf8("ysawtooth_radioButton"))
        self.gridLayout_2.addWidget(self.ysawtooth_radioButton, 0, 2, 1, 1)
        self.sawtoothMovement_lineEdit = QtGui.QLineEdit(self.groupBox_2)
        self.sawtoothMovement_lineEdit.setMaximumSize(QtCore.QSize(40, 16777215))
        self.sawtoothMovement_lineEdit.setObjectName(_fromUtf8("sawtoothMovement_lineEdit"))
        self.gridLayout_2.addWidget(self.sawtoothMovement_lineEdit, 0, 3, 1, 1)
        self.label = QtGui.QLabel(self.groupBox_2)
        self.label.setText(QtGui.QApplication.translate("Form", "microns", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 0, 4, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox_2)
        self.label_2.setText(QtGui.QApplication.translate("Form", "repeats", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_2.addWidget(self.label_2, 0, 6, 1, 1)
        self.sawtoothRepeats_lineEdit = QtGui.QLineEdit(self.groupBox_2)
        self.sawtoothRepeats_lineEdit.setMaximumSize(QtCore.QSize(40, 16777215))
        self.sawtoothRepeats_lineEdit.setObjectName(_fromUtf8("sawtoothRepeats_lineEdit"))
        self.gridLayout_2.addWidget(self.sawtoothRepeats_lineEdit, 0, 5, 1, 1)
        self.sawtoothTimeDelay_lineEdit = QtGui.QLineEdit(self.groupBox_2)
        self.sawtoothTimeDelay_lineEdit.setMaximumSize(QtCore.QSize(40, 16777215))
        self.sawtoothTimeDelay_lineEdit.setObjectName(_fromUtf8("sawtoothTimeDelay_lineEdit"))
        self.gridLayout_2.addWidget(self.sawtoothTimeDelay_lineEdit, 0, 7, 1, 1)
        self.label_3 = QtGui.QLabel(self.groupBox_2)
        self.label_3.setText(QtGui.QApplication.translate("Form", "ms time delay", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_2.addWidget(self.label_3, 0, 8, 1, 1)
        self.go_pushButton = QtGui.QPushButton(self.groupBox_2)
        self.go_pushButton.setText(QtGui.QApplication.translate("Form", "Go", None, QtGui.QApplication.UnicodeUTF8))
        self.go_pushButton.setObjectName(_fromUtf8("go_pushButton"))
        self.gridLayout_2.addWidget(self.go_pushButton, 1, 8, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_2)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

