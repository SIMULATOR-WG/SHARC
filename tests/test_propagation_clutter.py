# -*- coding: utf-8 -*-

"""
Created on Tue Mai 02 15:02:31 2017

@author: LeticiaValle_Mac
"""

import unittest
import numpy as np


from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss

class PropagationClutterLossTest(unittest.TestCase):

    def setUp(self):
        self.__ClutterAtt = PropagationClutterLoss(np.random.RandomState())

    def test_loss(self):

        f = 27000    #GHz
        theta = 0
        per_p =50
        d = 10000
        P = 0.9
        dist = 10
        r = 12.64
        s = 3.72
        t = 0.96
        u = 9.6
        v = 2
        w = 9.1
        x = -3
        y = 4.5
        z = -2
        #npt.assert_allclose(73.150,
         #                self.__ClutterAtt.get_loss(frequency=f, distance = d,percentage_p = per_p, dist = dist, elevation_angle_facade=theta, probability_loss_notExceeded=P, coeff_r=r, coeff_s=s, coeff_t=t, coeff_u=u,coeff_v=v, coeff_w=w,coeff_x=x,coeff_y=y,coeff_z=z),atol=1e-3)


#        f = [10,20]    #GHz
#        d = 10
#        Ph = 1013
#        T = 288
#        ro = 7.5
#        npt.assert_allclose([0.140, 1.088],
#                         self.__gasAtt.get_loss_Ag(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)


#        d = [[10, 20, 30],[40, 50, 60]]
#        f = 10
#        Ph = 1013
#        T = 288
#        ro = 7.5
#        self.assertTrue(np.all(np.isclose([0.140, 0.280, 0.420],[0.560, 0.700, 0.840],
#                        self.__gasAtt.get_loss_Ag(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro), atol=1e-3)))
#
#


if __name__ == '__main__':
    unittest.main()

