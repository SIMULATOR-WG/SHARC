# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:51:19 2017

@author: edgar
"""

from antenna import Antenna

class Station(object):
    
    def __init__(self):
        self.__id = -1
        self.__x = 0
        self.__y = 0
        self.__height = 0
        self.__tx_power = 0
        self.__rx_power = 0
        self.__tx_antenna = Antenna()
        self.__rx_antenna = Antenna()
        # test attributes below
        self.__active = False
        self.__acs = 0
        self.__aclr = 0
        self.__sensitivity = 0
        self.__noise_figure = 0
        self.__rx_aggr_interference = -300

    def set_id(self, id):
        self.__id = id

    def get_id(self):
        return self.__id

    def set_x(self, x):
        self.__x = x

    def get_x(self):
        return self.__x

    def set_y(self, y):
        self.__y = y

    def get_y(self):
        return self.__y

    def set_height(self, height):
        self.__height = height

    def get_height(self):
        return self.__height

    def set_tx_power(self, tx_power):
        self.__tx_power = tx_power

    def get_tx_power(self):
        return self.__tx_power

    def set_rx_power(self, rx_power):
        self.__rx_power = rx_power

    def get_rx_power(self):
        return self.__rx_power

    def set_tx_antenna(self, tx_antenna):
        self.__tx_antenna = tx_antenna

    def get_tx_antenna(self):
        return self.__tx_antenna

    def set_rx_antenna(self, rx_antenna):
        self.__rx_antenna = rx_antenna

    def get_rx_antenna(self):
        return self.__rx_antenna
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            equal = (self.get_id() == other.get_id() and 
                self.get_x() == other.get_x() and
                self.get_y() == other.get_y() and
                self.get_height() == other.get_height())
            return equal                
        else:
            return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            not_equal = (self.get_id() != other.get_id() or 
                self.get_x() != other.get_x() or
                self.get_y() != other.get_y() or
                self.get_height() != other.get_height())
            return not_equal                
        else:
            return NotImplemented  
            
    def __hash__(self):
        return hash((self.get_id(), self.get_x(), self.get_y(), 
                     self.get_height()))