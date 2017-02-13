# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:33:14 2017

@author: edgar
"""

from station import Station

class Link(object):
    
    def __init__(self):
        self.__transmitter = Station()
        self.__receiver = Station()
        self.__bandwidth = 0
        self.__frequency = 0
        self.__thermal_noise = 0
        self.__acir = 0
        self.__coupling_loss = 0