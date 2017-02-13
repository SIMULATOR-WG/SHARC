# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:38:20 2017

@author: edgar
"""

class Antenna(object):
    
    def __init__(self, gain=0):
        self.__gain = gain
        pass
    
    def set_gain(self, gain):
        self.__gain = gain
    
    def get_gain(self):
        return self.__gain

    def __float__(self):
        return float(self.get_gain())
        
    def __add__(self, other):
        return self.get_gain() + other

    def __radd__(self, other):
        return self.get_gain() + other

    def __sub__(self, other):
        return self.get_gain() - other

    def __rsub__(self, other):
        return other - self.get_gain()

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() < other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() <= other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() <= other
        else:
            return NotImplemented
        
    def __gt__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() > other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() > other
        else:
            return NotImplemented
    
    def __ge__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() >= other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() >= other
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() == other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() == other
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return self.get_gain() != other.get_gain()
        elif isinstance(other, (int, float)):
            return self.get_gain() != other            
        else:
            return NotImplemented            
