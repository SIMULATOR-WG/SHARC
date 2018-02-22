# -*- coding: utf-8 -*-
"""
Created on Tue Mai 02 15:02:31 2017

@author: LeticiaValle_Mac
"""

import unittest
import numpy as np

import numpy.testing as npt

from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

class PropagationGasesAttenuationTest(unittest.TestCase):

    def setUp(self):
        self.__gasAtt = PropagationGasesAttenuation(np.random.RandomState())

    def test_loss(self):

        f = 27    #GHz
        d = 10
        Ph = 1013
        T = 288
        ro = 7.5
        np.allclose(1.099,
                         self.__gasAtt.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)

        f = 27   #GHz
        d = [10,20]
        Ph = 1013
        T = 288
        ro = 7.5
        np.allclose([0.140, 3.297],
                         self.__gasAtt.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)


        f = 40   #GHz
        d = [10,20]
        Ph = 1013
        T = 288
        ro = 7.5
        np.allclose([1.295, 3.885],
                         self.__gasAtt.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)


if __name__ == '__main__':
    unittest.main()

