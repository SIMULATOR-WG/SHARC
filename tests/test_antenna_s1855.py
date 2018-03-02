#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 09:26:03 2017

@author: carlosrodriguez
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.antenna.antenna_s1855 import AntennaS1855
from sharc.parameters.parameters_fss_es import ParametersFssEs

class AntennaS1855Test(unittest.TestCase):

    def setUp(self):
        #Earth Station Antenna parameters
        params = ParametersFssEs()
        params.diameter = 9.1
        params.frequency = 27200
        params.antenna_gain = 62
        params.azimuth = 0
        params.elevation = 0

        # Create antenna FSS Earth Station objects
        self.antenna = AntennaS1855(params)

    def test_get_gain(self):

        off_axis_angle = np.array([7, 8, 15, 100])
        theta = np.array([90, 45, 45, 45])
        expected_result = np.array([ 10.87, 8.71, 2.59, -10 ])
        gain = self.antenna.calculate_gain(off_axis_angle_vec = off_axis_angle, 
                                           theta_vec = theta)
        npt.assert_allclose(gain, expected_result, atol=1e-2)

if __name__ == '__main__':
    unittest.main()
