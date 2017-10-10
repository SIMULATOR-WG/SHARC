# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 12:17:34 2017

@author: Andre Barreto
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.propagation.atmosphere import ReferenceAtmosphere

class TestAtmosphere(unittest.TestCase):

    def setUp(self):
        self.atmosphere = ReferenceAtmosphere()

    def test_specific_attenuation(self):
        temperature = 15 + 273.15 # K
        vapour_density = 7.5 # g/m**3
        pressure_hPa = 1013.25
        vapour_pressure_hPa = vapour_density * temperature / 216.7

        # compare with benchmark from ITU-R P-676-11 Figs. 1 and 2
        f_GHz_vec = [50, 60, 100, 200, 500, 1000]
        specific_att = np.zeros(len(f_GHz_vec))
        specific_att_p676_lower = [3e-1, 1e1, 4e-1, 2, 5e1, 6e2]
        specific_att_p676_upper = [5e-1, 2e1, 5e-1, 4, 7e1, 8e2]

        for index in range(len(f_GHz_vec)):
            specific_att[index] = self.atmosphere._get_specific_attenuation(pressure_hPa,
                                                                            vapour_pressure_hPa,
                                                                            temperature,
                                                                            float(f_GHz_vec[index]) * 1000)
        npt.assert_array_less(specific_att_p676_lower, specific_att)
        npt.assert_array_less(specific_att, specific_att_p676_upper)

if __name__ == '__main__':
    unittest.main()
