# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 12:17:34 2017

@author: Andre Barreto
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.propagation.propagation_p619 import PropagationP619


class DummyParams():
    def __init__(self):
        self.imt_altitude = 0
        self.surf_water_vapour_density = 0
        self.EARTH_RADIUS = 6371000


class TestPropagationP619(unittest.TestCase):

    def setUp(self):
        self.p619 = PropagationP619(np.random.RandomState())

    def test_atmospheric_gasses_loss (self):
        # compare with benchmark from ITU-R P-619 Fig. 3
        frequency_MHz = 30000.
        sat_params = DummyParams()
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
            loss = self.p619._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                          apparent_elevation=apparent_elevation,
                                                          surf_water_vapour_density=vapour_density,
                                                          sat_params=sat_params)
            self.assertLessEqual(loss_lower, loss)
            self.assertGreaterEqual(loss_upper, loss)


    def test_beam_spreading_attenuation(self):
        # compare with benchmark from ITU-R P-619 Fig. 7

        altitude_vec = np.array([0,1,2,3,4,6]) * 1000
        elevation_vec = [0,.5,1,2,3,5]
        att_lower_vec = [.8, .6 , .4, .2, .1, .05]
        att_upper_vec = [.9, .7, .5, .3, .2, .1]
        earth_to_space_vec = [True, False, True, False, True, False]

        for altitude, elevation, lower, upper, earth_to_space in zip(altitude_vec,
                                                                     elevation_vec,
                                                                     att_lower_vec,
                                                                     att_upper_vec,
                                                                     earth_to_space_vec):

            attenuation = self.p619._get_beam_spreading_att(elevation, altitude, earth_to_space)
            self.assertLessEqual(lower, abs(attenuation))
            self.assertGreaterEqual(upper, abs(attenuation))

            if earth_to_space:
                self.assertGreater(attenuation, 0)
            else:
                self.assertLess(attenuation, 0)


if __name__ == '__main__':
    unittest.main()
