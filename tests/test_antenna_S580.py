#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:44:40 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_s580 import AntennaS580
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np
import numpy.testing as npt

class AntennaS580Test(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssEs()
        param.antenna_gain = 40
        param.antenna_pattern = "ITU-R S580-6"
        param.diameter = 9.6
        param.frequency = 27250

        self.antenna = AntennaS580(param)
              
        
    def test_calculate_gain(self):
        psi = np.array([0,	1,	2,	3,	4,	5,	6,	7,	8,	10,	40, 60])

        ref_gain = np.array([40,	29,	21.47,	17.07,	13.94,	11.52,	9.54,	7.87,	6.42,	4,	-8.05,	-10])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    

        
if __name__ == '__main__':
    unittest.main()    