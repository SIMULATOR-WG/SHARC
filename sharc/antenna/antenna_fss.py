# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:18:59 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from parameters.parameters_fss import ParametersFss

class AntennaFss(Antenna):
    """
    Implements the satellite antenna pattern in the fixed-satellite service
    according to Recommendation ITU-R S.672-4
    """
    
    def __init__(self, param: ParametersFss):
        super(Antenna, self).__init__()
        self.peak_gain = param.sat_rx_antenna_gain
        
        
     