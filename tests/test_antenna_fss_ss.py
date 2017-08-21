# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 16:49:23 2017

@author: edgar
"""

import unittest

from sharc.antenna.antenna_fss_ss import AntennaFssSs
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaFssSsTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssSs()
        param.antenna_gain = 50
        param.antenna_pattern = "FSS_SS"
        param.antenna_3_dB = 2
        
        param.antenna_l_s = -25
        self.antenna25 = AntennaFssSs(param)
        
        param.antenna_l_s = -30    
        self.antenna30 = AntennaFssSs(param)
        
        
    def test_calculate_gain(self):
        psi = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100])

        ref_gain25 = np.array([0, -12, -28, -28, -28, -28, -28, -29.12, -30.57, -31.85, -33, -50])
        gain25 = self.antenna25.calculate_gain(phi_vec=psi) - self.antenna25.peak_gain
        npt.assert_allclose(gain25, ref_gain25, atol=1e-2)
    
        ref_gain30 = np.array([0, -12, -33, -33, -33, -33, -33, -34.12, -35.57, -36.85, -38, -50])
        gain30 = self.antenna30.calculate_gain(phi_vec=psi) - self.antenna30.peak_gain
        npt.assert_allclose(gain30, ref_gain30, atol=1e-2)

        
if __name__ == '__main__':
    unittest.main()    