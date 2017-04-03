# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:04:03 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
from support.observable import Observable
 
class Simulation(Observable):
    
    __metaclass__ = ABCMeta
    
    def __init__(self):
        super(Simulation, self).__init__()
    
    @abstractmethod
    def initialize(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def snapshot(self, *args, **kwargs):
        pass

    @abstractmethod
    def finalize(self, *args, **kwargs):
        pass
    