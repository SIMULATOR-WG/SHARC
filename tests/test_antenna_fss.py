# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:58:50 2017

@author: edgar
"""

import unittest

from sharc.antenna.antenna_fss import AntennaFss
from sharc.parameters.parameters_fss import ParametersFss

class AntennaFssTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFss()
        param.sat_rx_antenna_gain = 51
        param.sat_rx_antenna_pattern = "ITU-R S.672-4"
        param.sat_rx_antenna_l_s = -20    
        param.sat_rx_antenna_3_dB = 0.65
        
        self.antenna = AntennaFss(param)
        
    def test_gain(self):
        pass
    
        
if __name__ == '__main__':
    unittest.main()    