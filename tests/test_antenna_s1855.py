#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:26:03 2017

@author: carlosrodriguez
"""

import unittest
import numpy as np

from sharc.antenna.antenna_s1855 import AntennaS1855
#from sharc.parameters.parameters_antenna_fss import ParametersAntennaImt

class AntennaS1855Test(unittest.TestCase):
    
    def setUp(self):
        #Earth Station Antenna parameters
        self.diameter = 9.1
        self.frequency = 27200
        self.peak_gain = 62
        
        # Create antenna FSS Earth Station objects
        self.antenna = AntennaS1855(self.diameter, self.frequency, self.peak_gain)
                
    def test_get_gain(self):  
        # phi = 7 
        # phi = 8 with second part of equation
        # phi = 15, the third part of the equation
        # phi = 15, the third part of the equation
        phi = np.array([7, 8, 15, 100])
        theta = np.array([90, 45, 45, 45])
        expected_result = np.array([[10.872549, 9.372549, 9.372549, 9.372549], [9.53636364, 8.71818182, 8.71818182, 8.71818182], [2.5977185236079663, 2.5977185236079663, 2.5977185236079663, 2.5977185236079663], [-10.00, -10.00, -10.00, -10.00]])
        gain = self.antenna.get_gain(phi,theta)
        np.testing.assert_array_almost_equal(gain, expected_result)

        
if __name__ == '__main__':
    unittest.main()