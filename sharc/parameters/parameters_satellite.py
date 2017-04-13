# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:16:02 2017

@author: edgar
"""

class ParametersSatellite(object):

    __instance = None
    
    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersSatellite.__instance is None:
            ParametersSatellite.__instance = object.__new__(cls)
        return ParametersSatellite.__instance    
    
    ###########################################################################
    # 
    param = "VALUE"
