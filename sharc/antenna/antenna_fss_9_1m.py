# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:18:59 2017

@author: edgar
"""

import math
import numpy as np

from sharc.antenna.antenna import Antenna

class AntennaFss_9_1(Antenna):
    """"
    Implements amntenna radiation pattern for Earth Station according
    Model Viasat VA-91-KA used for 29 GHz 
    
    Attributes
    ----------
        gain (float): calculated antenna gain in given direction

    """
    
    def __init__(self,  peak_gain: float):
        """
        Constructs an AntennaFss_9_1 object.
        
        Parameters
        ---------
        peak_gain: peak gain of the antena on the direction of the GSO
        """

        self.__peak_gain = peak_gain
        
    def calculate_gain(self, *args, **kwargs) -> np.array: 
        """
        Calculates the gain of the antenna af arrays of angles.
            
        Returns
        -------
            gain (numpy.array): gain array in given directions
        """
        phi_list = np.arange (0,10,0.1)
        
        response = np.array([0,-4,-20,	-28,	-30,	-50,	-54,	-48,	-51,	-49,	-66,	-48,
                          -50,	-55,	-52,	-58,	-53,	-54,	-62,	-55,	-80,	-74,	-58,	-68,
                          -88,	-54,	-66,	-70,	-62,	-88,	-63,	-66,	-62,	-68,	-64,	-75,
                          -64,	-60,-90,	-60,	-82,	-69,	-64,	-67,	-90,	-72,	-90,	-72,	
                          -72,	-82,	-90,	-70,	-90,	-71,	-74,	-56,	-57,	-58,	-57,	-56,
                          -74,	-70,	-68,	-76,	-68,	-72,	-72,	-68,	-72,	-67,	-69,	-69,	
                          -69,	-69,	-70,	-72,	-74,	-76,	-78,	-88,	-76,	-75,	-90,	-90,	
                          -82,	-90,	-82,-62,	-86,	-74,	-90,	-78,	-88,	-90,	-88,	-90,
                          -87,	-90,	-88,	-82,	-88])
        return response       
     