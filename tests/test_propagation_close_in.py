# -*- coding: utf-8 -*-
"""
Created on Mon Mar  13 15:14:34 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.propagation.propagation_close_in import PropagationCloseIn

class PropagationCloseInTest(unittest.TestCase):

    def setUp(self):
        self.propagation = PropagationCloseIn(np.random.RandomState())

    def test_loss_los(self):
        d = np.array([10, 100])
        f = np.array([10, 100])
        los_prob = 1
        s_std = 0
        result = np.array([-27.55 + 20 + 20, -27.55 + 40 + 40])
        npt.assert_allclose(result, self.propagation.get_loss(distance_2D=d,
                                                              frequency=f,
                                                              line_of_sight_prob=los_prob,
                                                              shadowing=s_std))

        d = np.array([[1, 10],[100, 1000]])
        f = np.array([100])
        los_prob = 1
        s_std = 0
        result = np.array([[ -27.55 + 0 + 40,  -27.55 + 20 + 40 ], [ -27.55 + 40 + 40, -27.55 + 60 + 40]])
        npt.assert_allclose(result, self.propagation.get_loss(distance_2D=d,
                                                              frequency=f,
                                                              line_of_sight_prob=los_prob,
                                                              shadowing=s_std))


    def test_loss_nlos(self):
        d = np.array([10, 100])
        f = np.array([10, 100])
        los_prob = 0
        s_std = 0
        result = np.array([-27.55 + 30 + 20, -27.55 + 60 + 40])
        npt.assert_allclose(result, self.propagation.get_loss(distance_2D=d,
                                                              frequency=f,
                                                              line_of_sight_prob=los_prob,
                                                              shadowing=s_std))

        d = np.array([[1, 10],[100, 1000]])
        f = np.array([100])
        los_prob = 0
        s_std = 0
        result = np.array([[-27.55 + 0 + 40, -27.55 + 30 + 40], [-27.55 + 60 + 40, -27.55 + 90 + 40]])
        npt.assert_allclose(result, self.propagation.get_loss(distance_2D=d,
                                                              frequency=f,
                                                              line_of_sight_prob=los_prob,
                                                              shadowing=s_std))


    def test_shadowing(self):
        d = np.ones(100000)
        f = np.ones(100000)
        los_prob = 1
        s_std = 4.1
        loss = self.propagation.get_loss(distance_2D=d,
                                         frequency=f,
                                         line_of_sight_prob=los_prob,
                                         shadowing=s_std)

        mean = np.mean(loss)
        std = np.std(loss)

        self.assertAlmostEqual(mean, -27.55, delta = .1)
        self.assertAlmostEqual(std, 4.1, delta = .1)


if __name__ == '__main__':
    unittest.main()
