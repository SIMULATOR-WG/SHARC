# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:19:38 2017

@author: Calil
"""

import unittest

from sharc.antenna_imt import AntennaImt

class AntennaImtTest(unittest.TestCase):
    
    def setUp(self):
        #Element parameters
        g_max = 5
        theta_3db = 60
        phy_3db = 80
        am = 30
        sla_v = 30
        
        # Create antenna object
        self.antenna1 = AntennaImt(g_max,theta_3db,phy_3db,am,sla_v)
        
    def test_gain(self):
        self.assertEqual(self.antenna1.gain,0)
        
    def test_g_max(self):
        self.assertEqual(self.antenna1.g_max,5)
        
    def test_theta_3db(self):
        self.assertEqual(self.antenna1.theta_3db,60)
        
    def test_phy_3db(self):
        self.assertEqual(self.antenna1.phy_3db,80)
        
    def test_am(self):
        self.assertEqual(self.antenna1.am,30)
        
    def test_sla_v(self):
        self.assertEqual(self.antenna1.sla_v,30)
        
    def test_add(self):
        self.assertEqual(self.antenna1 + 2, 2)
    
    def test_radd(self):
        self.assertEqual(2 + self.antenna1, 2)
        
    def test_sub(self):
        self.assertEqual(self.antenna1 - 1, -1)
    
    def test_rsub(self):
        self.assertEqual(9 - self.antenna1, 9)
     
    def test_lt(self):
        self.assertTrue(self.antenna1 < 5)
        self.assertFalse(self.antenna1 < -3)
        
    def test_le(self):
        self.assertTrue(self.antenna1 <= 5)
        self.assertTrue(self.antenna1 <= 0)
        self.assertFalse(self.antenna1 <= -2)

    def test_gt(self):
        self.assertFalse(self.antenna1 > 5)
        self.assertTrue(self.antenna1 > -2)
        
    def test_ge(self):
        self.assertTrue(self.antenna1 >= 0)
        self.assertTrue(self.antenna1 >= -3)
        self.assertFalse(self.antenna1 >= 3)
        
    def test_eq(self):
        self.assertTrue(self.antenna1 == 0)
        self.assertFalse(self.antenna1 == 7)
        
    def test_ne(self):
        self.assertTrue(self.antenna1 != 1)
        self.assertFalse(self.antenna1 != 0)
        
    def test_vertical_pattern(self):
        # theta = 0 results in zero gain
        theta = 0
        v_att = self.antenna1.vertica_pattern(theta)
        self.assertEqual(v_att,0.0)
        
        # theta = 180 implies vertical gain of -27 dB
        theta = 180
        v_att = self.antenna1.vertica_pattern(theta)
        self.assertEqual(v_att,-27.0)
        
        # theta = 210, vertical attenuation equals vertical sidelobe attenuation
        theta = 210
        v_att = self.antenna1.vertica_pattern(theta)
        self.assertEqual(v_att,-30)
        self.assertEqual(v_att,-1.0*self.antenna1.sla_v)
        
    def test_horizontal_pattern(self):  
        # phy = 90 results in zero gain
        phy = 90
        h_att = self.antenna1.horizontal_pattern(phy)
        self.assertEqual(h_att,0.0)
        
        # phy = 120 implies horizontal gain of of -27 dB
        phy = 120
        h_att = self.antenna1.horizontal_pattern(phy)
        self.assertEqual(h_att,-27.0)
        
        # phy = 150, horizontal attenuation equals to the front-to-back ratio
        phy = 150
        h_att = self.antenna1.horizontal_pattern(phy)
        self.assertEqual(h_att,-30)
        self.assertEqual(h_att,-1.0*self.antenna1.am)
        
    def test_element_pattern(self):
        # theta = 0 and phy = 90 result in maximum gain
        theta = 0
        phy = 90
        e_gain = self.antenna1.element_pattern(theta,phy)
        self.assertEqual(e_gain,5.0)
        self.assertEqual(e_gain,self.antenna1.g_max)
        
        # Check for gain update
        self.assertTrue(self.antenna1 == 5.0)
        
        theta = 150
        phy = 80
        e_gain = self.antenna1.element_pattern(theta,phy)
        self.assertEqual(e_gain,-19.0)
        self.assertTrue(self.antenna1 == -19.0)
        
        theta = 210
        phy = 150
        e_gain = self.antenna1.element_pattern(theta,phy)
        self.assertEqual(e_gain,-25.0)
        self.assertEqual(e_gain,self.antenna1.g_max - self.antenna1.am)
        self.assertTrue(self.antenna1 == -25.0)
        
if __name__ == '__main__':
    unittest.main()