# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:56:47 2017

@author: edgar
"""

import unittest
import numpy.testing as npt
import numpy as np

from sharc.propagation.propagation_sat_simple import PropagationSatSimple


class PropagationSatSimpleTest(unittest.TestCase):

    def setUp(self):
        self.propagation = PropagationSatSimple(np.random.RandomState())

    def test_loss(self):
        d = np.array(10)
        f = np.array(10)
        loc_percentage = 0
        elevation = {'free_space': 90*np.ones(d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        self.assertEqual(12.45 + 4,
                         self.propagation.get_loss(distance_3D=d,
                                                   frequency=f,
                                                   loc_percentage=loc_percentage,
                                                   indoor_stations=indoor_stations,
                                                   elevation=elevation,
                                                   enable_clutter_loss=False))

        d = np.array([ 10, 100 ])
        f = np.array([ 10, 100 ])
        loc_percentage = 0
        elevation = {'free_space': 90*np.ones(d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        npt.assert_allclose([12.45 + 4, 52.45 + 4],
                            self.propagation.get_loss(distance_3D=d,
                                                      frequency=f,
                                                      loc_percentage=loc_percentage,
                                                      indoor_stations=indoor_stations,
                                                      elevation=elevation,
                                                      enable_clutter_loss=False),
                            atol=1e-2)

        d = np.array([ 10, 100, 1000 ])
        f = np.array([ 10, 100, 1000 ])
        loc_percentage = 0
        elevation = {'free_space': 90*np.ones(d.shape)}
        indoor_stations = np.array([1, 0, 0], dtype=bool)
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4, 92.45 + 4],
                            self.propagation.get_loss(distance_3D=d,
                                                      frequency=f,
                                                      loc_percentage=loc_percentage,
                                                      indoor_stations=indoor_stations,
                                                      elevation=elevation,
                                                      enable_clutter_loss=False),
                            atol=1e-2)

        d = np.array([[10, 20, 30],[40, 50, 60]])
        f = np.array([ 100 ])
        loc_percentage = 0
        elevation = {'free_space': 90*np.ones(d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        ref_loss = [[ 32.45 + 4,  38.47 + 4,  41.99 + 4],
                    [ 44.49 + 4,  46.42 + 4,  48.01 + 4]]
        loss = self.propagation.get_loss(distance_3D=d,
                                                   frequency=f,
                                                   loc_percentage=loc_percentage,
                                                   indoor_stations=indoor_stations,
                                                   elevation=elevation,
                                                   enable_clutter_loss=False)
        npt.assert_allclose(ref_loss, loss, atol=1e-2)


if __name__ == '__main__':
    unittest.main()
