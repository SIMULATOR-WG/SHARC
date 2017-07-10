# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:04:03 2017

@author: edgar
"""

from abc import ABC, abstractmethod
from sharc.support.observable import Observable
 
class Simulation(ABC, Observable):
    
    def __init__(self):
        ABC.__init__(self)
        Observable.__init__(self)
    
        
    @abstractmethod
    def initialize(self, *args, **kwargs):
        """
        This method is executed only once to initialize the simulation variables. 
        """        
        pass
    
    
    @abstractmethod
    def snapshot(self, *args, **kwargs):
        """
        Performs a single snapshot 
        """
        pass

    
    @abstractmethod
    def finalize(self, *args, **kwargs):
        """
        Finalizes the simulation (collect final results, etc...)
        """
        pass
    