# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:38:20 2017

@author: edgar
"""

class Antenna(object):
    """
    This is an omnidirectional antenna
    """
    
    def __init__(self, gain: float = 0):
        self.gain = gain
        
    """
    TODO: check the validity of operator overriding for Antenna class because
    it has to take into account the departure/arrival angle. Consider the 
    posibility of creating an arrival angle attribute (instead of parameter)
    """
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
