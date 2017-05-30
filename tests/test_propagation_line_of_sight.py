# -*- coding: utf-8 -*-


import unittest
import numpy as np

import numpy.testing as npt

from sharc.propagation.P452.propagation_line_of_sight import PropagationLineOfSight

class PropagationLineOfSightTest(unittest.TestCase):
    
    def setUp(self):
        self.__LineOfSight = PropagationLineOfSight()
        
    def test_loss(self):
        d = 10000
        f = 27000
        Ph = 1013 
        T = 288
        ro = 7.5
        
        npt.assert_allclose(142.226, 
                         self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)
        
        d = 20000
        f = 40000
        
        npt.assert_allclose(153.152, 
                         self.__LineOfSight.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)