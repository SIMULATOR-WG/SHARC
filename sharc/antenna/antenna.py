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
        self.beams_list = []
        self.w_vec_list = []
    
    
    @abstractmethod
    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the antenan gain.
        """        
        pass
    
    
    def add_beam(self, phi_etilt: float, theta_etilt: float):
        """
        Add new beam to antenna.
        Does not receive angles in local coordinate system.
        Theta taken with z axis as reference.
        
        Parameters
        ----------
            phi_etilt (float): azimuth electrical tilt angle [degrees]
            theta_etilt (float): elevation electrical tilt angle [degrees]
        """
        pass