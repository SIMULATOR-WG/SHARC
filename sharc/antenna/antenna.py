# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:09:17 2017

@author: Calil
"""

from abc import ABC, abstractmethod
import numpy as np

class Antenna(ABC):
    """
    Abstract antenna class. All antenna classes must inherit from it.
    
    Methods
    -------
    calculate_gain: calculates the antenna gain in the given directions
    """
    
    def __init__(self):
        pass
    
    
    @abstractmethod
    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the antenan gain.
        """        
        pass