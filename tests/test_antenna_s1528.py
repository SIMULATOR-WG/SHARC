# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 17:07:07 2017

@author: edgar
"""

import unittest

from sharc.antenna.antenna_s1528 import AntennaS1528
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaS1528Test(unittest.TestCase):

    def setUp(self):
        param = ParametersFssSs()
        param.antenna_gain = 39
        param.antenna_pattern = "ITU-R S.1528-0"
        param.antenna_3_dB = 2

        param.antenna_l_s = -20
        self.antenna20 = AntennaS1528(param)

        param.antenna_l_s = -30
        self.antenna30 = AntennaS1528(param)

    def test_calculate_gain(self):
        psi = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 80, 100])

        ref_gain20 = np.array([0, -3, -8.48, -20, -20, -20, -20, -21.10, -22.55, -23.83, -24.98, -39, -34.25])
        gain20 = self.antenna20.calculate_gain(off_axis_angle_vec=psi) - self.antenna20.peak_gain
        npt.assert_allclose(gain20, ref_gain20, atol=1e-2)

        ref_gain30 = np.array([0, -3, -8.48, -30, -30, -30, -30, -31.10, -32.55, -33.83, -34.98, -39, -39])
        gain30 = self.antenna30.calculate_gain(off_axis_angle_vec=psi) - self.antenna30.peak_gain
        npt.assert_allclose(gain30, ref_gain30, atol=1e-2)

if __name__ == '__main__':
    unittest.main()
