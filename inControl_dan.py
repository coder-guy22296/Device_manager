#!/usr/bin/python

#
#	A scientific instrumentation control software
#
#	Joerg 01/2012
#


import sys
sys.path.append('D:\\Programs')
import os
from PyQt4 import QtGui
import inLib
from inControl_ui import UI


class Control(object):
    '''
    The main class of InControl.
    It has each device and module as an attribute, as well as other instrument-wide properties.
    It also calls the main UI.
    '''
	
    def __init__(self):
        '''
        Reads the settings and loads all specified devices and modules.
        Sets up instrument-wide properties and calls the UI.
        '''
        print 'inControl: Starting up...'
        if len(sys.argv) > 1:
            settings_filename = sys.argv[1]
        else:
            print 'Loading default settings.'
            settings_filename = 'settings_default.yaml'
        self._settings = inLib.load_settings(settings_filename)

        # 1. Call device UI APIs

        #: A dictionary of all devices run by InControl.
        #: The keys are the name of a device, the value is an instance of the device's
        #: Control class. This dictionary is updated into the namespace of InControl's main
        #: control class. Devices are simple implementations to control one device.
        self._device_controls = {}
        if self._settings['devices']:
            for device in self._settings['devices']:
                device_dict = self._settings['devices'][device]
                if device_dict['active']:
                    device_path = inLib.get_device_path(device_dict)
                    device_module = inLib.import_module(device_path)
                    device_settings_dict = device_dict['settings']
                    device_settings_dict['settings_filename'] = settings_filename
                    device_control = device_module.Control(device_settings_dict)
                    self._device_controls[device] = inLib.CommandQueue(device_control)
        self.__dict__.update(self._device_controls)

        # 2. Call module UI APIs

        #: A dictionary of all modules run by InControl.
        #: The keys are the name of a module, the value is an instance of the module's
        #: Control class. This dictionary is updated into the namespace of InControl's main
        #: control class. InControl modules do not run a particular device, but have access to
        #: all device and other module control classes. Therefore, InControl's main control class
        #: is passed to the module object upon initiation. Modules should be used whenever a task
        #: with multiple devices needs to be performed.
        self._module_controls = {}
        for module in self._settings['modules']:
            if self._settings['modules'][module]['active']:
                module_path = inLib.get_module_path(module)
                module_module = inLib.import_module(module_path)
                module_settings = self._settings['modules'][module]['settings']
                module_control = module_module.Control(self, module_settings)
                self._module_controls[module] = inLib.CommandQueue(module_control)
        self.__dict__.update(self._module_controls)

        self.setWorkingDir(self._settings['working_dir'])

    
    def getDeviceControls(self):
        return self._device_controls


    def getMachineName(self):
        '''

        :Returns:
            *name*: str
        '''
        return self._settings['machine']


    def getModuleControls(self):
        return self._module_controls


    def getSettings(self):
        return self._settings


    def getWorkingDir(self):
        '''

        :Returns:
            *directory*: str
        '''
        return self._settings['working_dir']


    def setMachineName(self, name):
        '''
        Sets the name of the machine InControl is running on.

        :Parameters:
            *name*: str
        '''
        self._settings['machine'] = name


    def setWorkingDir(self, directory):
        '''
        Sets the current working directory of InControl and updates the respective label
        in the UI.

        :Parameters:
            *directory*: str
        '''
        self._settings['working_dir'] = directory
        os.chdir(directory)


    def shutDown(self):
        '''
        Shuts down all InControl modules and devices.
        It is called when the user closes InControl's main window.
        '''
        for module in self._module_controls.values():
            module.shutDown()
        for device in self._device_controls.values():
            device.shutDown()
        print 'Closing InControl.'

				
if __name__ == '__main__':
    control = Control()
    ui = UI(control)
