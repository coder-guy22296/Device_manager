# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'activationLED_design.ui'
#
# Created: Mon Feb 25 17:04:17 2013
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
        Form.resize(369, 420)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setAcceptDrops(True)
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Laser Control", None, QtGui.QApplication.UnicodeUTF8))
        self.horizontalLayout = QtGui.QHBoxLayout(Form)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.LED = QtGui.QGroupBox(Form)
        self.LED.setTitle(QtGui.QApplication.translate("Form", "LED", None, QtGui.QApplication.UnicodeUTF8))
        self.LED.setObjectName(_fromUtf8("LED"))
        self.gridLayout_5 = QtGui.QGridLayout(self.LED)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.OnLED = QtGui.QCheckBox(self.LED)
        font = QtGui.QFont()
        font.setPointSize(6)
        self.OnLED.setFont(font)
        self.OnLED.setText(QtGui.QApplication.translate("Form", "LED On", None, QtGui.QApplication.UnicodeUTF8))
        self.OnLED.setObjectName(_fromUtf8("OnLED"))
        self.gridLayout_5.addWidget(self.OnLED, 1, 0, 1, 1)
        self.SliderLED = QtGui.QSlider(self.LED)
        self.SliderLED.setMinimumSize(QtCore.QSize(0, 195))
        self.SliderLED.setMaximumSize(QtCore.QSize(16777215, 195))
        self.SliderLED.setOrientation(QtCore.Qt.Vertical)
        self.SliderLED.setObjectName(_fromUtf8("SliderLED"))
        self.gridLayout_5.addWidget(self.SliderLED, 0, 0, 1, 1)
        self.horizontalLayout.addWidget(self.LED)
        self.widget_shuttering = dropWidget(Form)
        self.widget_shuttering.setEnabled(True)
        self.widget_shuttering.setMouseTracking(False)
        self.widget_shuttering.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.widget_shuttering.setAcceptDrops(True)
        self.widget_shuttering.setObjectName(_fromUtf8("widget_shuttering"))
        self.formLayout = QtGui.QFormLayout(self.widget_shuttering)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label_11 = QtGui.QLabel(self.widget_shuttering)
        self.label_11.setText(QtGui.QApplication.translate("Form", "On time (sec):", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_11)
        self.lineEdit_onTime = QtGui.QLineEdit(self.widget_shuttering)
        self.lineEdit_onTime.setMaximumSize(QtCore.QSize(60, 16777215))
        self.lineEdit_onTime.setObjectName(_fromUtf8("lineEdit_onTime"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.lineEdit_onTime)
        self.label_12 = QtGui.QLabel(self.widget_shuttering)
        self.label_12.setText(QtGui.QApplication.translate("Form", "Period time (sec):", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_12)
        self.lineEdit_periodTime = QtGui.QLineEdit(self.widget_shuttering)
        self.lineEdit_periodTime.setMaximumSize(QtCore.QSize(60, 16777215))
        self.lineEdit_periodTime.setObjectName(_fromUtf8("lineEdit_periodTime"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.lineEdit_periodTime)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.formLayout.setItem(5, QtGui.QFormLayout.FieldRole, spacerItem)
        self.pushButton_start = QtGui.QPushButton(self.widget_shuttering)
        self.pushButton_start.setMaximumSize(QtCore.QSize(150, 16777215))
        self.pushButton_start.setText(QtGui.QApplication.translate("Form", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_start.setObjectName(_fromUtf8("pushButton_start"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.FieldRole, self.pushButton_start)
        self.label_cycles = QtGui.QLabel(self.widget_shuttering)
        self.label_cycles.setText(QtGui.QApplication.translate("Form", "Cycles:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_cycles.setObjectName(_fromUtf8("label_cycles"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_cycles)
        self.lineEdit_cycles = QtGui.QLineEdit(self.widget_shuttering)
        self.lineEdit_cycles.setMaximumSize(QtCore.QSize(60, 16777215))
        self.lineEdit_cycles.setObjectName(_fromUtf8("lineEdit_cycles"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.lineEdit_cycles)
        self.horizontalLayout.addWidget(self.widget_shuttering)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        pass

from myWidget.dragAndDrop import dropWidget
