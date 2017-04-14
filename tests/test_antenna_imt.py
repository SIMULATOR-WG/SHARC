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
        theta_3db = 65
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
        self.assertEqual(self.antenna1.theta_3db,65)
        
    def test_phy_3db(self):
        self.assertEqual(self.antenna1.phy_3db,80)
        
    def test_am(self):
        self.assertEqual(self.antenna1.am,30)
        
    def test_sla_v(self):
        self.assertEqual(self.antenna1.sla_v,30)
        
if __name__ == '__main__':
    unittest.main()