# -*- coding: utf-8 -*-

"""
Created on  Mar  04 14:13:31 2017

@author:  LeticiaValle_Mac
"""

import unittest
import numpy as np
import numpy.testing as npt


from sharc.propagation.propagation_abg import PropagationABG

class PropagationABGTest(unittest.TestCase):

    def setUp(self):
        self.abg = PropagationABG(np.random.RandomState())

    def test_loss(self):
        d = np.array([[100, 500],[400, 60]])
        f = 27000
        indoor = np.zeros(d.shape, dtype=bool)
        alpha = 3.4
        beta = 19.2
        gamma = 2.3
        shadowing = 0
        loss = np.array ([[120.121, 143.886347],[140.591406, 112.578509]])

        npt.assert_allclose(self.abg.get_loss(distance_2D = d,
                                              frequency = f,
                                              indoor_stations = indoor,
                                              line_of_sight_prob = 1,
                                              alpha = alpha,
                                              beta = beta,
                                              gamma = gamma,
                                              shadowing = shadowing),
                             loss,atol=1e-2)

        d = np.array([500, 3000])
        f = np.array([27000, 40000])
        indoor = np.zeros(d.shape, dtype=bool)
        alpha = 3.4
        beta = 19.2
        gamma = 2.3
        shadowing = 0

        loss = np.array ([143.886,174.269])
        npt.assert_allclose(self.abg.get_loss(distance_2D = d,
                                              frequency = f,
                                              indoor_stations = indoor,
                                              line_of_sight_prob = 1,
                                              alpha = alpha,
                                              beta = beta,
                                              gamma = gamma,
                                              shadowing = shadowing),
                           loss ,atol=1e-2)
