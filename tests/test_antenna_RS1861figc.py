#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:52:58 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_rs1861_fig9c import AntennaRS1861FIG9c
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaRS1861FIG9cTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssSs()
        param.antenna_pattern = "ITU-R RS.1869-figure c"
        param.antenna_gain = 50


        self.antenna = AntennaRS1861FIG9c(param)
        
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1, 2, 3, 4, 5, 6, 7, 8, 9, 	10,	100, 150])

        ref_gain = np.array([50,	48.42,	46.75,	44.97,	43.07,	40.99,	38.71,	36.11,	33.04,	29.07,	17,	-5.84, -16])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    
        
if __name__ == '__main__':
    unittest.main()   