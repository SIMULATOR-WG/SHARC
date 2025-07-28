# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:39:34 2017

@author: edgar
"""

import unittest
import numpy.testing as npt
import math

from sharc.antenna.antenna_rra7_3 import AntennaReg_RR_A7_3
from sharc.parameters.parameters_antenna_with_diameter import ParametersAntennaWithDiameter


class AntennaReg_RR_A7_3Test(unittest.TestCase):

    def setUp(self):
        self.antenna1_params = ParametersAntennaWithDiameter(
                diameter=3,
                frequency=10000,
                antenna_gain=50,
        )

        self.antenna1 = AntennaReg_RR_A7_3(
            self.antenna1_params,
        )
        self.antenna2_params = ParametersAntennaWithDiameter(
                diameter=3,
                frequency=4000,
                antenna_gain=30,
        )
        self.antenna2 = AntennaReg_RR_A7_3(
            self.antenna2_params,
        )
        self.antenna3_params = ParametersAntennaWithDiameter(
            frequency=10000,
            antenna_gain=50,
        )
        self.antenna3 = AntennaReg_RR_A7_3(
            self.antenna3_params,
        )

    def test_diameter_inference(self):
        """Document specifies a diameter to assume for when no diameter is given"""
        self.assertAlmostEqual(
            20 * math.log10(self.antenna3.D_lmbda),
            self.antenna3_params.antenna_gain - 7.7,
            delta=1e-10,
        )

    def test_gain(self):
        self.assertEqual(self.antenna1.peak_gain, self.antenna1_params.antenna_gain)
        self.assertEqual(self.antenna1.D_lmbda, 100)
        self.assertEqual(self.antenna1.g1, 29)
        self.assertAlmostEqual(self.antenna1.phi_r, 1.000067391, delta=1e-11)
        self.assertAlmostEqual(self.antenna1.phi_m, 0.91651513899, delta=1e-11)

        self.assertEqual(self.antenna2.peak_gain, self.antenna2_params.antenna_gain)
        self.assertEqual(self.antenna2.D_lmbda, 40)
        self.assertAlmostEqual(self.antenna2.g1, 19.0514997832, delta=1e-11)
        self.assertAlmostEqual(self.antenna2.phi_r, 2.5, delta=1e-11)
        self.assertAlmostEqual(self.antenna2.phi_m, 1.65442589867, delta=1e-11)

    def test_invalid_antenna(self):
        """
        The appendix 7 Annex 3 does not specify what should happen
        in case phi_r < phi_m. We're throwing an error, but if you
        happen to know that isn't an error, remove it and verify the behavior you want
        """

        with self.assertRaises(ValueError):
            AntennaReg_RR_A7_3(
                ParametersAntennaWithDiameter(
                    diameter=3,
                    frequency=10000,
                    antenna_gain=100,
                ),
            )

    def test_calculate_gain(self):
        # Test antenna1
        phi = [
            0.9,
            0.92,
            1.1,
            35.9,
            36,
            180,
        ]
        expected = [
            self.antenna1.peak_gain - 2.5 * 1e-3 * self.antenna1.D_lmbda * phi[0] * self.antenna1.D_lmbda * phi[0],
            self.antenna1.g1,
            29 - 25 * math.log10(phi[2]),
            29 - 25 * math.log10(phi[3]),
            -10,
            -10,
        ]

        # Test antenna2
        phi = [
            1.6,
            1.7,
            2.4,
            2.6,
            35.9,
            36,
            180,
        ]
        expected = [
            self.antenna2.peak_gain - 2.5 * 1e-3 * self.antenna2.D_lmbda * phi[0] * self.antenna2.D_lmbda * phi[0],
            self.antenna2.g1,
            self.antenna2.g1,
            29 - 25 * math.log10(phi[3]),
            29 - 25 * math.log10(phi[4]),
            -10,
            -10,
        ]

        # Test antenna1
        gains2 = self.antenna2.calculate_gain(off_axis_angle_vec=phi)
        self.assertEqual(len(gains2), len(phi))
        npt.assert_allclose(gains2, expected)


if __name__ == '__main__':
    unittest.main()
