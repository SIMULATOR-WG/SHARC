# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:38:20 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
import numpy as np

class AntennaOmni(Antenna):
    """
    This is an omnidirectional antenna
    """
    
    def __init__(self, gain: float = 0):
        self.gain = gain
    
    def calculate_gain(self,directions: list) -> np.array:
        """
        Calculates the gain, which is the same for all the directions
        """
        return self.gain*np.ones(len(directions))
    
    def __float__(self):
        return float(self.gain)
        
    def __add__(self, other):
        return self.gain + other

    def __radd__(self, other):
        return self.gain + other

    def __sub__(self, other):
        return self.gain - other

    def __rsub__(self, other):
        return other - self.gain

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.gain < other.gain
        elif isinstance(other, (int, float)):
            return self.gain < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, self.__class__):
            return self.gain <= other.gain
        elif isinstance(other, (int, float)):
            return self.gain <= other
        else:
            return NotImplemented
        
    def __gt__(self, other):
        if isinstance(other, self.__class__):
            return self.gain > other.gain
        elif isinstance(other, (int, float)):
            return self.gain > other
        else:
            return NotImplemented
    
    def __ge__(self, other):
        if isinstance(other, self.__class__):
            return self.gain >= other.gain
        elif isinstance(other, (int, float)):
            return self.gain >= other
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.gain == other.gain
        elif isinstance(other, (int, float)):
            return self.gain == other
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return self.gain != other.gain
        elif isinstance(other, (int, float)):
            return self.gain != other            
        else:
            return NotImplemented            
