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
        self.__gasAtt = PropagationGasesAttenuation()
        
    def test_loss(self):
       
        f = 10    #GHz
        d = 10
        Ph = 1013 
        T = 288
        ro = 7.5
        npt.assert_allclose(0.140, 
                         self.__gasAtt.get_loss_Ag(distance=d, frequency=f, atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro),atol=1e-3)
        
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
       
