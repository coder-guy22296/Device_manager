#!/usr/bin/python

#
#	A scientific instrumentation control software
#
#	Joerg 01/2012
#


import sys
from PyQt4 import QtGui
import inLib
import inControl_design
	
class UI(object):

    def __init__(self, control):
        self._control = control
        
        self._app = QtGui.QApplication(sys.argv)

        self._window = QtGui.QWidget()
        self._window.closeEvent = self.shutDown

        self._ui = inControl_design.Ui_Form()
        self._ui.setupUi(self._window)

        self._window.setWindowTitle('InControl: ' + self._control.getMachineName())
        self._ui.labelDir.setText(self._control.getWorkingDir())
        self._ui.lineEdit_workingDir.setText(self._control.getWorkingDir())

        self._ui.pushButton_setDir.clicked.connect(self.setWorkingDirFromLineEdit)

        self._ui.listWidgetDevices.itemActivated.connect(self._on_row_activatedD)
        self._ui.listWidgetModules.itemActivated.connect(self._on_row_activatedM)
        self._ui.pushButtonDir.clicked.connect(self.setWorkingDir)

        settings = self._control.getSettings()

        device_controls = self._control.getDeviceControls()

        self._device_ui_controls = {}

        # Check if settings contain a device section:
        if settings['devices']:
            # Loop over all devices:
            for device in settings['devices']:
                # Check if device is marked as active:
                if settings['devices'][device]['active']:
                    # Get the path to a potential ui module:
                    device_ui_path = inLib.get_device_ui_path(settings['devices'][device])
                    # Try to import the ui module:
                    try:
                        print 'Trying to import', device_ui_path
                        device_ui_module = inLib.import_module(device_ui_path)
                        #print "ui_path:", device_ui_path
                        # Start the UI:
                        device_control = device_controls[device]
                        self._device_ui_controls[device] = device_ui_module.UI(device_control)
                    except ImportError:
                        # There is no ui for this device. Do nothing.
                        print 'Did not find UI for {0}. Passing...'.format(device)
                    except:
                        raise
                    # Add an entry to the device list widget:
                    self._ui.listWidgetDevices.addItem(device)
        self.__dict__.update(self._device_ui_controls)

        module_controls = self._control.getModuleControls()
        self._module_ui_controls = {}
        for module in settings['modules']:
            if settings['modules'][module]['active']:
                module_ui_path = inLib.get_module_ui_path(module)
                try:
                    print 'Trying to import', module_ui_path
                    module_ui_module = inLib.import_module(module_ui_path)
                    self._module_ui_controls[module] = module_ui_module.UI(module_controls[module],
                            self)
                except ImportError:
                    print 'Did not find UI for {0}. Passing...'.format(module)
                self._ui.listWidgetModules.addItem(module)
        self.__dict__.update(self._module_ui_controls)

        self._window.show()
        self._app.exec_()


    def setWorkingDir(self):
        wd = self._control.getWorkingDir()
        working_dir = str(QtGui.QFileDialog.getExistingDirectory(self._window,
                                    'Select directory',
                                    wd))
        if working_dir:
            self._ui.labelDir.setText(working_dir)
            self._control.setWorkingDir(working_dir)

    def setWorkingDirFromLineEdit(self):
        wd = str(self._ui.lineEdit_workingDir.text())
        self._control.setWorkingDir(wd)
        self._ui.labelDir.setText(wd)


    def _on_row_activatedD(self, item):
        temp = str(item.text())
        if temp in self._device_ui_controls.keys():
            self._device_ui_controls[temp]._window.show()


    def _on_row_activatedM(self, item):
        temp = str(item.text())
        if temp in self._module_ui_controls.keys():
            self._module_ui_controls[temp]._window.show()


    def shutDown(self, event):
        for m in self._module_ui_controls.values():
            m.shutDown()
        for d in self._device_ui_controls.values():
            d.shutDown()
        self._control.shutDown()
        self._app.quit()
