#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:26:03 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_s1855 import AntennaS1855
#from sharc.parameters.parameters_antenna_fss import ParametersAntennaImt

class AntennaS1855Test(unittest.TestCase):
    
    def setUp(self):
        #Earth Station Antenna parameters
        self.diameter = 9.1
        self.frequency = 27200
        
        # Create antenna FSS Earth Station objects
        self.antenna = AntennaS1855(self.diameter, self.frequency)
                
    def test_get_gain(self):  
        # phi = 7 
        phi = 7
        theta = 90
        expected_result = 10.872549
        gain = self.antenna.get_gain(phi,theta)
        self.assertAlmostEqual(gain, expected_result, places=2)
        
        
        # phi = 8 with second part of equation
        phi = 8
        theta = 45
        expected_result = 8.718181818181819
        gain = self.antenna.get_gain(phi,theta)
        self.assertAlmostEqual(gain, expected_result, places=2)
        
        # phi = 15, the third part of the equation
        phi = 15
        expected_result = 2.5977185236079663
        gain = self.antenna.get_gain(phi,theta)
        self.assertAlmostEqual(gain, expected_result, places =2)

        # phi = 15, the third part of the equation
        phi = 100
        gain = self.antenna.get_gain(phi,theta)
        self.assertAlmostEqual(gain, -10.00, places=2)

        
if __name__ == '__main__':
    unittest.main()