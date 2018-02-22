
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 11:14:56 2017

@author: LeticiaValle_Mac
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.propagation.propagation_umi import PropagationUMi

class PropagationUMiTest(unittest.TestCase):

    def setUp(self):
        self.umi = PropagationUMi(np.random.RandomState())

    def test_los_probability(self):
        distance_2D = np.array([[10, 15, 40],
                                [17, 60, 80]])
        los_probability = np.array([[1, 1,  0.631],
                                    [1, 0.432, 0.308]])
        npt.assert_allclose(self.umi.get_los_probability(distance_2D), 
                            los_probability,
                            atol=1e-2)


    def test_breakpoint_distance(self):
        h_bs = np.array([15, 20, 25, 30])
        h_ue = np.array([3, 4])
        h_e = np.ones((len(h_bs), len(h_ue)))
        frequency = 30000*np.ones((len(h_bs), len(h_ue)))
        breakpoint_distance = np.array([[ 11200,  16800],
                                        [ 15200,  22800],
                                        [ 19200,  28800],
                                        [ 23200,  34800]])
        npt.assert_array_equal(self.umi.get_breakpoint_distance(frequency, h_bs, h_ue, h_e),
                               breakpoint_distance)


    def test_loss_los(self):
        distance_2D = np.array([[100, 200, 300, 400],
                                [500, 600, 700, 800]])
        h_bs = np.array([30, 35])
        h_ue = np.array([2, 3, 4, 5])
        h_e = np.ones(distance_2D.shape)
        distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
        frequency = 30000*np.ones(distance_2D.shape)
        shadowing_std = 0
        loss = np.array([[104.336, 110.396, 114.046, 116.653],
                         [118.690,  120.346,  121.748,  122.963]])
        npt.assert_allclose(self.umi.get_loss_los(distance_2D, distance_3D, frequency,
                                                  h_bs, h_ue, h_e, shadowing_std),
                            loss,
                            atol=1e-2)


        distance_2D = np.array([[100, 200, 300, 400],
                                [500, 600, 700, 800]])
        h_bs = np.array([30, 35])
        h_ue = np.array([2, 3, 4, 5])
        h_e = np.ones(distance_2D.shape)
        distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
        frequency = 300*np.ones(distance_2D.shape)
        shadowing_std = 0
        loss = np.array([[64.336, 70.396, 74.046, 76.653],
                         [89.215, 86.829, 86.187, 86.139]])
        npt.assert_allclose(self.umi.get_loss_los(distance_2D, distance_3D, frequency,
                                                  h_bs, h_ue, h_e, shadowing_std),
                            loss,
                            atol=1e-2)


    def test_loss_nlos(self):
        distance_2D = np.array([[100, 200, 300, 400],
                                [500, 600, 700, 800]])
        h_bs = np.array([30, 35])
        h_ue = np.array([2, 3, 4, 5])
        h_e = np.ones(distance_2D.shape)
        distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
        frequency = 30000*np.ones(distance_2D.shape)
        shadowing_std = 0
        loss = np.array([[128.841, 138.727, 144.562, 148.645],
                         [152.969, 155.453, 157.509, 159.252]])
        npt.assert_allclose(self.umi.get_loss_nlos(distance_2D, distance_3D, frequency,
                                                  h_bs, h_ue, h_e, shadowing_std),
                            loss,
                            atol=1e-2)


        distance_2D = np.array([[1000, 2000, 5000, 4000],
                                [3000, 6000, 7000, 8000]])
        h_bs = np.array([30, 35])
        h_ue = np.array([2, 3, 4, 5])
        h_e = np.ones(distance_2D.shape)
        distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
        frequency = 300*np.ones(distance_2D.shape)
        shadowing_std = 0
        loss = np.array([[120.968, 131.290, 145.036, 141.315],
                         [137.805, 148.131, 150.194, 151.941]])
        npt.assert_allclose(self.umi.get_loss_nlos(distance_2D, distance_3D, frequency,
                                                  h_bs, h_ue, h_e, shadowing_std),
                            loss,
                            atol=1e-2)


if __name__ == '__main__':
    unittest.main()
