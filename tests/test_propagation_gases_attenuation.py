# -*- coding: utf-8 -*-
"""
Created on Tue Mai 02 15:02:31 2017

@author: LeticiaValle_Mac
"""

import unittest
import numpy as np

import numpy.testing as npt

from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

class PropagationFreeSpaceTest(unittest.TestCase):
    
    def setUp(self):
        self.__gasAtt = PropagationGasesAttenuation()
        
    def test_loss(self):
       
        d = 10
        f = 10000    #MHz
        Ph = 1013 
        T = 288
        ro = 7.5
        npt.assert_allclose(72.590, 
                         self.__gasAtt.get_loss(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)


        d = [[10, 20, 30],[40, 50, 60]]
        f = [ 10000 ]
        Ph = 1013 
        T = 288
        ro = 7.5
        self.assertTrue(np.all(np.isclose([ 72.590,  78.610,  82.132],
                    [ 84.631,  86.569,  88.153], self.__gasAtt.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro), atol=1e-2)))   
   
    
    
    
if __name__ == '__main__':
    unittest.main()
       
