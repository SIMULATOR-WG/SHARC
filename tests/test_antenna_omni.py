# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:39:34 2017

@author: edgar
"""

import unittest
import numpy.testing as npt

from sharc.antenna.antenna_omni import AntennaOmni

class AntennaOmniTest(unittest.TestCase):
    
    def setUp(self):
        self.antenna1 = AntennaOmni()
        self.antenna1.gain = 5
        
        self.antenna2 = AntennaOmni()
        self.antenna2.gain = 8.0
        
        self.antenna3 = AntennaOmni(10)
        
    def test_gain(self):
        self.assertEqual(self.antenna1.gain, 5)
        self.assertEqual(self.antenna2.gain, 8)
        self.assertEqual(self.antenna3.gain, 10)
        
    def test_calculate_gain(self):
        phi = [30, 60, 90, 45]
        
        # Test antenna1
        gains1 = self.antenna1.calculate_gain(phi_vec=phi)
        self.assertEqual(len(gains1), len(phi))
        npt.assert_allclose(gains1, self.antenna1.gain)
        
        # Test antenna2
        gains2 = self.antenna2.calculate_gain(phi_vec=phi)
        self.assertEqual(len(gains2), len(phi))
        npt.assert_allclose(gains2, self.antenna2.gain)
        
        # Test antenna3
        gains3 = self.antenna3.calculate_gain(phi_vec=phi)
        self.assertEqual(len(gains3), len(phi))
        npt.assert_allclose(gains3, self.antenna3.gain)
        
        
if __name__ == '__main__':
    unittest.main()
