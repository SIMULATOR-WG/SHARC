# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 18:32:30 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.simulation_uplink import SimulationUplink
from sharc.parameters.parameters_imt import ParametersImt
from sharc.antenna.antenna_omni import AntennaOmni

class SimulationUplinkTest(unittest.TestCase):
    
    def setUp(self):
        self.param = ParametersImt()
        self.param.topology = "SINGLE_BS"
        self.param.num_base_stations = 1
        self.param.num_clusters = 1
        self.param.static_base_stations = True
        self.param.intersite_distance = 2000
        self.param.interfered_with = False
        self.param.frequency = 27500
        self.param.bandwidth = 100    
        self.param.mcl = 117
        self.param.ho_margin = 3
        self.param.bs_load_probability = 0.5
        self.param.bs_tx_power = 40
        self.param.bs_height = 10
        self.param.bs_tx_antenna_gain = 0
        self.param.bs_rx_antenna_gain = 0
        self.param.bs_aclr = 40    
        self.param.bs_acs = 30    
        self.param.bs_noise_figure = 7  
        self.param.bs_feed_loss = 3
        self.param.ue_k = 1
        self.param.ue_k_m = 1
        self.param.ue_tx_power_control = "OFF"
        self.param.ue_tx_power = 22
        self.param.ue_height = 1.5
        self.param.ue_tx_antenna_gain = 0
        self.param.ue_rx_antenna_gain = 0
        self.param.ue_aclr = 35    
        self.param.ue_acs = 25    
        self.param.ue_noise_figure = 9    
        self.param.ue_feed_loss = 3
        self.param.channel_model = "FSPL"
        
    def test_simulation_2bs_4ue(self):
        self.param.num_base_stations = 1
        self.param.num_clusters = 2
        self.param.ue_k = 2
        self.param.ue_k_m = 1
        self.simulation_uplink = SimulationUplink(self.param)
        
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_uplink.ue), 4)
        self.assertEqual(len(self.simulation_uplink.bs), 2)
        self.assertEqual(self.simulation_uplink.coupling_loss.shape, (2,4))
        
        # after initialize(), receivers (base stations) must be created
        self.simulation_uplink.initialize()
        self.assertEqual(self.simulation_uplink.bs.num_stations, 2)
        
        # create FSS system
        self.simulation_uplink.create_system()
        
        # it is time to create user equipments
        self.simulation_uplink.create_ue()
        self.simulation_uplink.ue.x = np.array([-2000, -500, 400, 1500])
        self.simulation_uplink.ue.y = np.array([0, 0, 0, 0])
        self.assertEqual(self.simulation_uplink.ue.num_stations, 4)

        self.simulation_uplink.bs.rx_antenna = [AntennaOmni(0), AntennaOmni(1)]
        self.simulation_uplink.ue.tx_antenna = [AntennaOmni(2), AntennaOmni(3), AntennaOmni(4), AntennaOmni(5)]
        
        # let's calculate coupling loss
        self.simulation_uplink.coupling_loss =  np.transpose( \
            self.simulation_uplink.calculate_coupling_loss(self.simulation_uplink.ue,
                                                           self.simulation_uplink.bs,
                                                           self.simulation_uplink.propagation_imt) )
        npt.assert_allclose(self.simulation_uplink.coupling_loss, 
                            np.array([[119.23,  112.21,  120.15,  124.19], [127.77,  120.75,  111.80,  109.21]]), 
                            atol=1e-2)
        
        # Now we connect base stations to user equipments
        self.simulation_uplink.connect_ue_to_bs()
        self.assertEqual(self.simulation_uplink.link, {0: [0,1], 1: [2,3]})
        
        self.simulation_uplink.select_ue()
        
        # Scheduling algorirhm
        self.simulation_uplink.scheduler()
        npt.assert_equal(self.simulation_uplink.ue.bandwidth,
                         45*np.ones(4))
        
        # apply power control to set transmit powers
        self.simulation_uplink.power_control()
        npt.assert_allclose(self.simulation_uplink.ue.tx_power, 
                            22*np.ones(4), 
                            atol=1e-2)

#        self.simulation_uplink.calculate_external_interference()   
#        self.assertAlmostEqual(self.simulation_uplink.system.inr, 1.02, places=2)
        

    def test_calculate_gains(self):
        self.param.num_base_stations = 1
        self.param.num_clusters = 2
        self.param.ue_k = 2
        self.param.ue_k_m = 1
        self.simulation_uplink = SimulationUplink(self.param)
        
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_uplink.ue), 4)
        self.assertEqual(len(self.simulation_uplink.bs), 2)
        self.assertEqual(self.simulation_uplink.coupling_loss.shape, (2,4))
        
        # after initialize(), receivers (base stations) must be created
        self.simulation_uplink.initialize()
        self.assertEqual(self.simulation_uplink.bs.num_stations, 2)
        
        # it is time to create user equipments
        self.simulation_uplink.create_ue()
        self.simulation_uplink.ue.x = np.array([-2000, -500, 400, 1500])
        self.simulation_uplink.ue.y = np.array([0, 0, 0, 0])
        self.assertEqual(self.simulation_uplink.ue.num_stations, 4)

        self.simulation_uplink.bs.rx_antenna = [AntennaOmni(0), AntennaOmni(1)]
        self.simulation_uplink.ue.tx_antenna = [AntennaOmni(2), AntennaOmni(3),\
                                                AntennaOmni(4), AntennaOmni(5)]
        
        # Now we calculate the gain matrix
        gains = self.simulation_uplink.calculate_gains(self.simulation_uplink.ue,\
                                                       self.simulation_uplink.bs,\
                                                       "TX")
        npt.assert_equal(gains,np.array([[2, 2],
                                         [3, 3],
                                         [4, 4],
                                         [5, 5]]))
        gains = self.simulation_uplink.calculate_gains(self.simulation_uplink.bs,\
                                                       self.simulation_uplink.ue,\
                                                       "RX")
        npt.assert_equal(gains,np.array([[0, 0, 0, 0],
                                         [1, 1, 1, 1]]))
                
if __name__ == '__main__':
    unittest.main()
          