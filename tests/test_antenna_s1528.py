# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 17:07:07 2017

@author: edgar
"""

import unittest

from sharc.antenna.antenna_s1528 import AntennaS1528, AntennaS1528Taylor
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528

import numpy as np
import numpy.testing as npt


class AntennaS1528Test(unittest.TestCase):

    def setUp(self):
        param = ParametersAntennaS1528()
        param.antenna_gain = 39
        param.antenna_pattern = "ITU-R S.1528-0"
        param.antenna_3_dB_bw = 2

        param.antenna_l_s = -20
        self.antenna20 = AntennaS1528(param)

        param.antenna_l_s = -30
        self.antenna30 = AntennaS1528(param)

    def test_calculate_gain(self):
        psi = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 80, 100])

        ref_gain20 = np.array([
            0, -
            3, -
            8.48, -
            20, -
            20, -
            20, -
            20, -
            21.10, -
            22.55, -
            23.83, -
            24.98, -
            39, -
            34.25,
        ])
        gain20 = self.antenna20.calculate_gain(
            off_axis_angle_vec=psi,
        ) - self.antenna20.peak_gain
        npt.assert_allclose(gain20, ref_gain20, atol=1e-2)

        ref_gain30 = np.array(
            [0, -3, -8.48, -30, -30, -30, -30, -31.10, -32.55, -33.83, -34.98, -39, -39],
        )
        gain30 = self.antenna30.calculate_gain(
            off_axis_angle_vec=psi,
        ) - self.antenna30.peak_gain
        npt.assert_allclose(gain30, ref_gain30, atol=1e-2)

    def test_calculate_params_bessel(self):
        """Compare the parameteres calculated by the class with the reference values present in the Recommendation
        S.1528-5 - Annex II - Examples for recommends 1.4
        """
        params_rolloff_7 = ParametersAntennaS1528(
            antenna_gain=0,
            frequency=12000,
            bandwidth=0,
            slr=20,
            n_side_lobes=4,
        )

        # Create an instance of AntennaS1528Taylor
        antenna_rolloff_7 = AntennaS1528Taylor(params_rolloff_7)
        ref_primary_roots = np.array([1.219, 2.233, 3.238])
        ref_A = 0.95277
        ref_sigma = 1.1692
        npt.assert_allclose(antenna_rolloff_7.mu, ref_primary_roots, atol=1e-3)
        npt.assert_allclose(antenna_rolloff_7.A, ref_A, atol=1e-5)
        npt.assert_allclose(antenna_rolloff_7.sigma, ref_sigma, atol=1e-4)


if __name__ == '__main__':
    unittest.main()
