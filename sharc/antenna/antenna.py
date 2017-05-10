# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:09:17 2017

@author: Calil
"""

from abc import ABCMeta, abstractmethod

class Antenna(object):
    """
    Abstract antenna class. All antenna classes must inherit from it.
    
    Methods
    -------
    calculate_gain: calculates the antenna gain in the 
        given directions
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def calculate_gain(self, *args, **kwargs):
        pass