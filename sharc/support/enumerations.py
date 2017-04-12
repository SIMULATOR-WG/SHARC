# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 17:53:41 2017

@author: edgar
"""

from enum import Enum

class Action( Enum ):
    """
    The action that is sent to controller in order to control the simulation
    """
    START_SIMULATION = 1
    STOP_SIMULATION  = 2
    
class State( Enum ):
    """
    This is the graphical user interface state
    """
    INITIAL  = 1
    RUNNING  = 2
    FINISHED = 3
    STOPPED  = 4
    STOPPING = 5