# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:09:17 2017

@author: Calil
"""

from abc import ABCMeta, abstractmethod
import numpy as np

class Antenna(object):
    """
    Abstract antenna class. All antenna classes must inherit from it.
    
    Methods
    -------
    calculate_gain(directions: np.array): calculates the antenna gain in the 
        given directions
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def calculate_gain(self,directions: np.array) -> np.array:
        """
        Calculates the gain in the given direction.
        
        Parameters
        ----------
        directions (np.array): array of tuples, each containing the aximuth 
            (phi) and elevation (theta) angles to which the gain is calculated.
            
        Returns
        -------
        gains (np.array): gain corresponding to each of the given directions.
        """
        pass