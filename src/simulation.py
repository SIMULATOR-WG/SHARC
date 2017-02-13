# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:04:03 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
 
class Simulation(object):
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def initialize(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def snapshot(self, *args, **kwargs):
        pass

    @abstractmethod
    def finalize(self, *args, **kwargs):
        pass
    