#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:26:03 2017

@author: carlosrodriguez
"""

import unittest
import numpy as np

from sharc.antenna.antenna_fss_9_1m import AntennaFss_9_1
#from sharc.parameters.parameters_antenna_fss import ParametersAntennaImt

class AntennaFSS_9_1mTest(unittest.TestCase):
    
    def setUp(self):
        #Earth Station Antenna parameters
        self.diameter = 9.1
        self.frequency = 27200
        self.peak_gain = 62
        
        # Create antenna FSS Earth Station objects
        self.antenna = AntennaFss_9_1(self.peak_gain)
                
    def test_get_gain(self):  

        phi = np.arange (0,10,0.1)
        expected_result = np.array([0,-4,-20,	-28,	-30,	-50,	-54,	-48,	-51,	-49,	-66,	-48,
                          -50,	-55,	-52,	-58,	-53,	-54,	-62,	-55,	-80,	-74,	-58,	-68,
                          -88,	-54,	-66,	-70,	-62,	-88,	-63,	-66,	-62,	-68,	-64,	-75,
                          -64,	-60,-90,	-60,	-82,	-69,	-64,	-67,	-90,	-72,	-90,	-72,	
                          -72,	-82,	-90,	-70,	-90,	-71,	-74,	-56,	-57,	-58,	-57,	-56,
                          -74,	-70,	-68,	-76,	-68,	-72,	-72,	-68,	-72,	-67,	-69,	-69,	
                          -69,	-69,	-70,	-72,	-74,	-76,	-78,	-88,	-76,	-75,	-90,	-90,	
                          -82,	-90,	-82,-62,	-86,	-74,	-90,	-78,	-88,	-90,	-88,	-90,
                          -87,	-90,	-88,	-82,	-88])
        #gain = self.antenna.get_gain()
        #np.testing.assert_array_almost_equal(gain, expected_result)

        
if __name__ == '__main__':
    unittest.main()