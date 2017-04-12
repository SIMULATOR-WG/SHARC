# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:05:46 2017

@author: edgar
"""

class ParametersGeneral(object):
        
    __instance = None
    
    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """        
        if ParametersGeneral.__instance is None:
            ParametersGeneral.__instance = object.__new__(cls)
        return ParametersGeneral.__instance    
    
    ###########################################################################
    # Number of simulation snapshots
    num_snapshots = 100
