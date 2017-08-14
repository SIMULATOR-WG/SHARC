# -*- coding: utf-8 -*-
"""
Created on Mon Mar  13 15:14:34 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.parameters.parameters_fss import ParametersFss
import matplotlib.pyplot as plt

class TestPropagationP619(unittest.TestCase):

    def setUp(self):
        self.p619 = PropagationP619()

    def test_specific_attenuation(self, plotFlag = False):
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
            specific_att[index] = self.p619._get_specific_attenuation(pressure_hPa,
                                                                      vapour_pressure_hPa,
                                                                      temperature,
                                                                      float(f_GHz_vec[index]) * 1000)
        npt.assert_array_less(specific_att_p676_lower, specific_att)
        npt.assert_array_less(specific_att, specific_att_p676_upper)

        if plotFlag:
            # generate plot
            f_GHz_vec = range(1,1000)
            specific_att = np.zeros(len(f_GHz_vec))

            for index in range(len(f_GHz_vec)):
                specific_att[index] = self.p619._get_specific_attenuation(pressure_hPa,
                                                                          vapour_pressure_hPa,
                                                                          temperature,
                                                                          float(f_GHz_vec[index]) * 1000)
            plt.figure()
            plt.semilogy( f_GHz_vec, specific_att )
            plt.xlabel('frequency(GHz)')
            plt.xlabel('Specific attenuation (dB/km)')

    def test_atmospheric_gasses_loss (self, plotFlag = False):
        # compare with benchmark from ITU-R P-619 Fig. 3
        frequency_MHz = 30000.
        sat_params = ParametersFss()
        sat_params.imt_altitude = 1000

        vapour_density_vec = [7.5, 12.5, 2.5, 5.]
        apparent_elevation_vec = [0, 10, 20, 40]

        loss_lower_vec = [10, 1, .3, .2]
        loss_upper_vec = [20, 2, .4, .3]

        for vapour_density, apparent_elevation, loss_lower, loss_upper in zip(vapour_density_vec,
                                                                              apparent_elevation_vec,
                                                                              loss_lower_vec,
                                                                              loss_upper_vec):
            sat_params.surf_water_vapour_density = vapour_density
            loss = self.p619._get_atmospheric_gasses_loss(frequency_MHz, apparent_elevation,
                                                          sat_params)
            self.assertLessEqual(loss_lower, loss)
            self.assertGreaterEqual(loss_upper, loss)

        if plotFlag:
            apparent_elevation = range(-1, 90, 2)
            loss_2_5 = np.zeros(len(apparent_elevation))
            loss_12_5 = np.zeros(len(apparent_elevation))

            for index in range(len(apparent_elevation)):
                print("Apparent Elevation: {} degrees".format(apparent_elevation[index]))

                sat_params.surf_water_vapour_density = 2.5
                loss_2_5[index] = self.p619._get_atmospheric_gasses_loss(frequency_MHz,
                                                                         apparent_elevation[index],
                                                                         sat_params)
                sat_params.surf_water_vapour_density = 12.5
                loss_12_5[index] = self.p619._get_atmospheric_gasses_loss(frequency_MHz,
                                                                         apparent_elevation[index],
                                                                         sat_params)

            plt.figure()
            plt.semilogy(apparent_elevation, loss_2_5, label='2.5 g/m^3')
            plt.semilogy(apparent_elevation, loss_12_5, label='12.5 g/m^3')

            plt.xlabel("apparent elevation (deg)")
            plt.ylabel("Loss (dB)")
            plt.legend()


if __name__ == '__main__':
    unittest.main()
