""" 
This is the control program for coherent obis 488 laser. 

Created by: Dan Xie, 07/27/2016.
"""
import inLib
import sys 
import os 
import time


class Control(inLib.Device):
    def __init__(self, settings):
        self.executable = settings['executable']
        self.status = False
        
        
    def turn_on(self):
        pass 
    
    def stand_by(self):
        pass

