#!/usr/bin/python

import inLib
import sys
import traceback
import imp


class Control(inLib.Module):

    def __init__(self, control, settings):
        print 'Initializing Scripts.'
        inLib.Module.__init__(self, control, settings)
        self.stop = False

    def run(self, script):
        '''
        Runs a script and passes an instance of the InControls main class, so that devices
        and modules can be access in the script by *control*.

        :Parameters:
            *script*: str
                The sscript.
        '''
        if False:
            try:
                exec(script, {'control' : self._control})
            except:
                etype, value, tb = sys.exc_info()
                print ''.join(traceback.format_exception(etype, value, tb))
        else:
            # Indent the script:
            lines = script.splitlines(True)
            script = 'def run(control):\n'
            for line in lines:
                script = script + '    ' + line
            # Create a python module and execute the function:
            module = imp.new_module('_script_module')
            exec script in module.__dict__
            # Run the function and pass inControl' control instance:
            module.run(self._control)
