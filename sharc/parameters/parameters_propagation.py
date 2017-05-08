# -*- coding: utf-8 -*-
"""
Created on Fri May  5 16:36:47 2017

@author: edgar
"""


class ParametersPropagation(object):

    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersPropagation.__instance is None:
            ParametersPropagation.__instance = object.__new__(cls)
        return ParametersPropagation.__instance


    ###########################################################################
    # Total air pressure in hPa
    atmospheric_pressure = 1013
    
    ###########################################################################
    # Temperature in Kelvin
    air_temperature = 288
    
    ###########################################################################
    # water vapour concentration (g/m^3)
    water_vapour = 7.5