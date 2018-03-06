# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 10:30:15 2018

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.antenna.beamforming_normalization.beamforming_normalization import BeamformingNormalization

class BeamformingNormalizationTest(unittest.TestCase):
    
    def setUp(self):
        # Test 1
        resolution = 30
        self.norm_1 = BeamformingNormalization(resolution)
    
    def test_construction(self):
        # Test 1
        self.assertEqual(self.norm_1.res,30)
        self.assertEqual(self.norm_1.phi_min,-180)
        self.assertEqual(self.norm_1.phi_max,180)
        self.assertEqual(self.norm_1.theta_min,0)
        self.assertEqual(self.norm_1.theta_max,180)
        npt.assert_equal(self.norm_1.phi_vals,np.arange(-180,180,30))
        npt.assert_equal(self.norm_1.theta_vals,np.arange(0,180,30))
    
    def test_generate_correction_matrix(self):
        pass
    
    def test_calculate_correction_factor(self):
        pass
    
    def test_save_files(self):
        pass
    
if __name__ == '__main__':
    unittest.main()