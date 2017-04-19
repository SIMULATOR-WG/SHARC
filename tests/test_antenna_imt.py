# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:19:38 2017

@author: Calil
"""

import unittest

from sharc.antenna.antenna_imt import AntennaImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaImtTest(unittest.TestCase):
    
    def setUp(self):
        #Element parameters
        self.param = ParametersAntennaImt()
        self.param.bs_element_max_g = 5
        self.param.bs_element_phi_3db = 80
        self.param.bs_element_theta_3db = 60
        self.param.bs_element_am = 30
        self.param.bs_element_sla_v = 30
        self.param.ue_element_max_g = 10
        self.param.ue_element_phi_3db = 75
        self.param.ue_element_theta_3db = 65
        self.param.ue_element_am = 25
        self.param.ue_element_sla_v = 35
        
        # Create antenna IMT objects
        self.antenna1 = AntennaImt(self.param,"BS")
        self.antenna2 = AntennaImt(self.param,"UE")
        self.antenna2.gain = 15
        
    def test_gain(self):
        self.assertEqual(self.antenna1.gain,0)
        self.assertEqual(self.antenna2.gain,15)
        
    def test_station_type(self):
        self.assertTrue(self.antenna1.station_type == "BS")
        self.assertTrue(self.antenna2.station_type == "UE")
        
    def test_g_max(self):
        self.assertEqual(self.antenna1.g_max,5)
        self.assertEqual(self.antenna2.g_max,10)
        
    def test_phi_3db(self):
        self.assertEqual(self.antenna1.phi_3db,80)
        self.assertEqual(self.antenna2.phi_3db,75)
        
    def test_theta_3db(self):
        self.assertEqual(self.antenna1.theta_3db,60)
        self.assertEqual(self.antenna2.theta_3db,65)
        
    def test_am(self):
        self.assertEqual(self.antenna1.am,30)
        self.assertEqual(self.antenna2.am,25)
        
    def test_sla_v(self):
        self.assertEqual(self.antenna1.sla_v,30)
        self.assertEqual(self.antenna2.sla_v,35)
        
    def test_add(self):
        self.assertEqual(self.antenna1 + 2, 2)
    
    def test_radd(self):
        self.assertEqual(2 + self.antenna1, 2)
        self.assertEqual(self.antenna1 + self.antenna2, 15)
        
    def test_sub(self):
        self.assertEqual(self.antenna1 - 1, -1)
    
    def test_rsub(self):
        self.assertEqual(9 - self.antenna1, 9)
        self.assertEqual(self.antenna2 - self.antenna1, 15)
     
    def test_lt(self):
        self.assertTrue(self.antenna1 < 5)
        self.assertFalse(self.antenna1 < -3)
        self.assertTrue(self.antenna1 < self.antenna2)
        
    def test_le(self):
        self.assertTrue(self.antenna1 <= 5)
        self.assertTrue(self.antenna1 <= 0)
        self.assertFalse(self.antenna1 <= -2)
        self.assertTrue(self.antenna1 <= self.antenna2)

    def test_gt(self):
        self.assertFalse(self.antenna1 > 5)
        self.assertTrue(self.antenna1 > -2)
        self.assertFalse(self.antenna1 > self.antenna2)
        
    def test_ge(self):
        self.assertTrue(self.antenna1 >= 0)
        self.assertTrue(self.antenna1 >= -3)
        self.assertFalse(self.antenna1 >= 3)
        self.assertFalse(self.antenna1 >= self.antenna2)
        
    def test_eq(self):
        self.assertTrue(self.antenna1 == 0)
        self.assertFalse(self.antenna1 == 7)
        self.assertFalse(self.antenna1 == self.antenna2)
        
    def test_ne(self):
        self.assertTrue(self.antenna1 != 1)
        self.assertFalse(self.antenna1 != 0)
        self.assertTrue(self.antenna1 != self.antenna2)
        
    def test_horizontal_pattern(self):  
        # phi = 0 results in zero gain
        phi = 0
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,0.0)
        
        # phi = 120 implies horizontal gain of of -27 dB
        phi = 120
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,-27.0)
        
        # phi = 150, horizontal attenuation equals to the front-to-back ratio
        phi = 150
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,-30)
        self.assertEqual(h_att,-1.0*self.antenna1.am)
        
    def test_vertical_pattern(self):
        # theta = 90 results in zero gain
        theta = 90
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,0.0)
        
        # theta = 180 implies vertical gain of -27 dB
        theta = 180
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,-27.0)
        
        # theta = 210, vertical attenuation equals vertical sidelobe attenuation
        theta = 210
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,-30)
        self.assertEqual(v_att,-1.0*self.antenna1.sla_v)
        
    def test_element_pattern(self):
        # theta = 0 and phi = 90 result in maximum gain
        phi = 0
        theta = 90
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,5.0)
        self.assertEqual(e_gain,self.antenna1.g_max)
        
        # Check for gain update
        self.assertTrue(self.antenna1 == 5.0)
        
        phi = 80
        theta = 150
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-19.0)
        self.assertTrue(self.antenna1 == -19.0)
        
        phi = 150
        theta = 210
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-25.0)
        self.assertEqual(e_gain,self.antenna1.g_max - self.antenna1.am)
        self.assertTrue(self.antenna1 == -25.0)
        
if __name__ == '__main__':
    unittest.main()