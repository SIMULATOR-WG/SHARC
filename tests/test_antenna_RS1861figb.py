#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:47:43 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_rs1861_fig9b import AntennaRS1861FIG9b
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaRS1861FIG9bTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssSs()
        param.antenna_pattern = "ITU-R RS.1869-figure b"
        param.antenna_gain = 50


        self.antenna = AntennaRS1861FIG9b(param)
        
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1, 2, 3, 4, 5, 6, 7, 8, 9, 	10,	100])

        ref_gain = np.array([50,	29.99,	11.35,	8.05,	4.76,	1.46,	-1.82,	-5.11,	-8.41,	-11.70,	-15,	-15])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    
        
if __name__ == '__main__':
    unittest.main()   