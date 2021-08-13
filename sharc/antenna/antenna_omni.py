# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:38:20 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
import numpy as np

class AntennaOmni(Antenna):
    """
    This is an omnidirectional antenna
    """

    def __init__(self, gain: float = 0):
        super().__init__()
        self.gain = gain

    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the gain, which is the same for all the directions

        Parameters
        ----------
        phi_vec (np.array): azimuth angles [degrees]

        Returns
        -------
        gains (np.array): numpy array of gains
        """

        if "phi_vec" in kwargs:
            phi_vec = np.asarray(kwargs["phi_vec"])
        else:
            phi_vec = np.asarray(kwargs["off_axis_angle_vec"])

        return self.gain*np.ones(len(phi_vec))

