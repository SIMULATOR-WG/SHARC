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
        self.propagation = PropagationSatSimple(
            random_number_gen=np.random.RandomState(),
            enable_clutter_loss=False)

    def test_loss(self):
        d = np.array(10)
        f = np.array(10)
        elevation = {
            'free_space': 90 *
            np.ones(
                d.shape),
            'apparent': 90 *
            np.ones(
                d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        self.assertAlmostEqual(13.19,
                               self.propagation.get_loss(
                                   d, f, indoor_stations, elevation),
                               delta=1e-2)

        d = np.array([10, 100])
        f = np.array([10, 100])
        elevation = {
            'free_space': 90 *
            np.ones(
                d.shape),
            'apparent': 90 *
            np.ones(
                d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        npt.assert_allclose([13.2, 53.2],
                            self.propagation.get_loss(
                                d, f, indoor_stations, elevation),
                            atol=1e-2)

        d = np.array([10, 100, 1000])
        f = np.array([10, 100, 1000])
        elevation = {
            'free_space': 90 *
            np.ones(
                d.shape),
            'apparent': 90 *
            np.ones(
                d.shape)}
        indoor_stations = np.array([0, 0, 0], dtype=bool)
        npt.assert_allclose([13.2, 53.2, 93.2],
                            self.propagation.get_loss(
                                d, f, indoor_stations, elevation),
                            atol=1e-2)

        d = np.array([[10, 20, 30], [40, 50, 60]])
        f = np.array([100])
        elevation = {
            'free_space': 90 *
            np.ones(
                d.shape),
            'apparent': 90 *
            np.ones(
                d.shape)}
        indoor_stations = np.zeros(d.shape, dtype=bool)
        ref_loss = [[33.2, 39.2, 42.7],
                    [45.2, 47.1, 48.7]]
        loss = self.propagation.get_loss(d, f, indoor_stations, elevation)
        npt.assert_allclose(ref_loss, loss, atol=1e-1)


if __name__ == '__main__':
    unittest.main()
