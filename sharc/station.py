# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:51:19 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna

class Station(object):
    
    def __init__(self):
        self.__id = -1
        self.__x = 0
        self.__y = 0
        self.__height = 0
        self.__tx_power = 0
        self.__rx_power = 0
        self.__antenna = Antenna()
        # test attributes below
        self.__active = False
        self.__acs = 0
        self.__aclr = 0
        self.__sensitivity = 0
        self.__noise_figure = 0
        self.__rx_aggr_interference = -300

    @property
    def id(self):
        return self.__id
        
    @id.setter
    def id(self, value):
        self.__id = value
        
    @property
    def x(self):
        return self.__x
        
    @x.setter
    def x(self, value):
        self.__x = value

    @property
    def y(self):
        return self.__y
        
    @y.setter
    def y(self, value):
        self.__y = value
        
    @property
    def height(self):
        return self.__height
        
    @height.setter
    def height(self, value):
        self.__height = value
        
    @property
    def active(self):
        return self.__active
        
    @active.setter
    def active(self, value):
        self.__active = value
        
    @property
    def tx_power(self):
        return self.__tx_power
        
    @tx_power.setter
    def tx_power(self, value):
        self.__tx_power = value

    @property
    def rx_power(self):
        return self.__rx_power
        
    @rx_power.setter
    def rx_power(self, value):
        self.__rx_power = value

    @property
    def antenna(self):
        return self.__antenna
        
    @antenna.setter
    def antenna(self, value):
        self.__antenna = value
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            equal = (self.id == other.id and 
                self.x == other.x and
                self.y == other.y and
                self.height == other.height)
            return equal
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            not_equal = (self.id != other.id or 
                self.x != other.x or
                self.y != other.y or
                self.height != other.height)
            return not_equal
        else:
            return NotImplemented
            
    def __hash__(self):
        return hash(self.id, self.x, self.y, self.height)
