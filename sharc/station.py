# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:51:19 2017

@author: edgar
"""

from sharc.support.enumerations import StationType

class Station(object):
    
    def __init__(self):
        self.id = -1
        self.x = 0
        self.y = 0
        self.azimuth = 0
        self.elevation = 0
        self.height = 0
        self.indoor = False
        self.active = False
        self.tx_power = 0
        self.rx_power = 0
        self.rx_interference = 0
        self.ext_interference = 0
        self.antenna = None
        self.bandwidth = 0
        self.noise_figure = 0
        self.noise_temperature = 0
        self.thermal_noise = 0
        self.total_interference = 0
        self.snr = 0
        self.sinr = 0
        self.sinr_ext = 0
        self.inr = 0
        self.station_type = StationType.NONE
        
        
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
