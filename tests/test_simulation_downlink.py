# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 16:22:42 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters_imt import ParametersImt
from sharc.antenna.antenna import Antenna

class SimulationDownlinkTest(unittest.TestCase):
    
    def setUp(self):
        self.param = ParametersImt()
        self.param.topology = "SINGLE_BS"
        self.param.num_base_stations = 1
        self.param.num_clusters = 1
        self.param.static_base_stations = True
        self.param.intersite_distance = 2000
        self.param.interfered_with = True
        self.param.frequency = 50000
        self.param.bandwidth = 200    
        self.param.mcl = 126
        self.param.ho_margin = 3
        self.param.bs_load_probability = 0.5
        self.param.num_resource_blocks = 10
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
        self.param.ue_tx_power = 22
        self.param.ue_height = 1.5
        self.param.ue_tx_antenna_gain = 0
        self.param.ue_rx_antenna_gain = 0
        self.param.ue_aclr = 35    
        self.param.ue_acs = 25    
        self.param.ue_noise_figure = 9    
        self.param.ue_feed_loss = 3
        
    def test_simulation_1bs_1ue(self):
        self.simulation_downlink = SimulationDownlink(self.param)
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_downlink.transmitter), 1)
        self.assertEqual(len(self.simulation_downlink.receiver), 1)
        self.assertEqual(self.simulation_downlink.coupling_loss.shape, (1,1))
        
        # after initialize(), transmitters must be created
        self.simulation_downlink.initialize()
        self.assertEqual(self.simulation_downlink.transmitter.num_stations, 1)
        
        # it is time to create user equipments
        self.simulation_downlink.create_ue()
        self.simulation_downlink.receiver.x = np.array([100])
        self.simulation_downlink.receiver.y = np.array([0])
        self.assertEqual(self.simulation_downlink.receiver.num_stations, 1)
        
        # let's calculate coupling loss
        self.simulation_downlink.calculate_coupling_loss()
        npt.assert_allclose(self.simulation_downlink.coupling_loss, np.array([[106.43]]), atol=1e-2)
        
        self.simulation_downlink.connect_ue_to_bs()
        self.assertEqual(self.simulation_downlink.link, {0: [0]})
        
        self.simulation_downlink.scheduler()
        npt.assert_equal(self.simulation_downlink.receiver.bandwidth, 200)
        
        self.simulation_downlink.power_control()
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[0], 
                            [40], atol=1e-2)
        
        self.simulation_downlink.calculate_sinr()      
        npt.assert_allclose(self.simulation_downlink.receiver.rx_power,
                            [-66.43], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.rx_interference,
                            [-300], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.thermal_noise,
                            -111.96, atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.total_interference,
                            [-111.96], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.sinr,
                            [45.54], atol=1e-2)
        
    def test_simulation_1bs_2ue(self):
        self.param.ue_k = 2
        self.simulation_downlink = SimulationDownlink(self.param)
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_downlink.transmitter), 1)
        self.assertEqual(len(self.simulation_downlink.receiver), 2)
        self.assertEqual(self.simulation_downlink.coupling_loss.shape, (1,2))
        
        # after initialize(), transmitters must be created
        self.simulation_downlink.initialize()
        self.assertEqual(self.simulation_downlink.transmitter.num_stations, 1)
        
        # it is time to create user equipments
        self.simulation_downlink.create_ue()
        self.simulation_downlink.receiver.x = np.array([100, 200])
        self.simulation_downlink.receiver.y = np.array([0, 0])
        self.assertEqual(self.simulation_downlink.receiver.num_stations, 2)
        
        # let's calculate coupling loss
        self.simulation_downlink.calculate_coupling_loss()
        npt.assert_allclose(self.simulation_downlink.coupling_loss, np.array([[106.43, 112.45]]), atol=1e-2)
        
        self.simulation_downlink.connect_ue_to_bs()
        self.assertEqual(self.simulation_downlink.link, {0: [0,1]})
        
        self.simulation_downlink.scheduler()
        npt.assert_equal(self.simulation_downlink.receiver.bandwidth, 100*np.ones(2))
        
        self.simulation_downlink.power_control()
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[0], 
                            [36.98,36.98], atol=1e-2)
        
        self.simulation_downlink.calculate_sinr()      
        npt.assert_allclose(self.simulation_downlink.receiver.rx_power,
                            [-69.43, -75.46], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.rx_interference,
                            [-300, -300], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.thermal_noise,
                            -114.97, atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.total_interference,
                            [-114.97, -114.97], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.sinr,
                            [45.54, 39.52], atol=1e-2)
                

    def test_simulation_1bs_4ue(self):
        self.param.ue_k = 4
        self.param.ue_k_m = 1
        self.simulation_downlink = SimulationDownlink(self.param)
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_downlink.transmitter), 1)
        self.assertEqual(len(self.simulation_downlink.receiver), 4)
        self.assertEqual(self.simulation_downlink.coupling_loss.shape, (1,4))
        
        # after initialize(), transmitters must be created
        self.simulation_downlink.initialize()
        self.assertEqual(self.simulation_downlink.transmitter.num_stations, 1)
        
        # it is time to create user equipments
        self.simulation_downlink.create_ue()
        self.simulation_downlink.receiver.x = np.array([100, 200, 400, 800])
        self.simulation_downlink.receiver.y = np.array([0, 0, 0, 0])
        self.assertEqual(self.simulation_downlink.receiver.num_stations, 4)
        
        # let's calculate coupling loss
        self.simulation_downlink.calculate_coupling_loss()
        npt.assert_allclose(self.simulation_downlink.coupling_loss, np.array([[106.43, 112.45, 118.47, 124.49]]), atol=1e-2)
        
        # Now we connect base stations to user equipments
        self.simulation_downlink.connect_ue_to_bs()
        self.assertEqual(self.simulation_downlink.link, {0: [0,1,2,3]})
        
        # Scheduling algorirhm
        self.simulation_downlink.scheduler()
        npt.assert_equal(self.simulation_downlink.receiver.bandwidth,
                         50*np.ones(4))
        
        # apply power control to set transmit powers
        self.simulation_downlink.power_control()
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[0], 
                            [33.97, 33.97, 33.97, 33.97], 
                            atol=1e-2)
        
        # calculate SINR
        self.simulation_downlink.calculate_sinr()
        npt.assert_allclose(self.simulation_downlink.receiver.rx_power,
                            [-72.45, -78.47, -84.49, -90.51],
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.rx_interference,
                            [-300,  -300,  -300, -300], 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.thermal_noise,
                            -117.98*np.ones(4), 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.total_interference,
                            -117.98*np.ones(4), 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.sinr,
                            [45.53,  39.51,  33.49,  27.47], 
                            atol=1e-2)
        
    def test_simulation_2bs_2ue(self):
        self.param.num_base_stations = 1
        self.param.num_clusters = 2
        self.param.mcl = 127
        
        self.simulation_downlink = SimulationDownlink(self.param)
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_downlink.transmitter), 2)
        self.assertEqual(len(self.simulation_downlink.receiver), 2)
        self.assertEqual(self.simulation_downlink.coupling_loss.shape, (2,2))
        
        # after initialize(), transmitters must be created
        self.simulation_downlink.initialize()
        self.assertEqual(self.simulation_downlink.transmitter.num_stations, 2)
        
        # it is time to create user equipments
        self.simulation_downlink.create_ue()
        self.simulation_downlink.receiver.x = np.array([-1200, 1700])
        self.simulation_downlink.receiver.y = np.array([0, 0])
        self.assertEqual(self.simulation_downlink.receiver.num_stations, 2)

        self.simulation_downlink.transmitter.tx_antenna = [Antenna(0), Antenna(1)]
        self.simulation_downlink.receiver.rx_antenna = [Antenna(2), Antenna(3)]
        
        # let's calculate coupling loss
        self.simulation_downlink.calculate_coupling_loss()
        npt.assert_allclose(self.simulation_downlink.coupling_loss, np.array([[110.45, 132.05],[130.27, 119.33]]), atol=1e-2)
        
        self.simulation_downlink.connect_ue_to_bs()
        self.assertEqual(self.simulation_downlink.link, {0: [0], 1:[1]})
        
        self.simulation_downlink.scheduler()
        npt.assert_equal(self.simulation_downlink.receiver.bandwidth, 200*np.ones(2))
        
        self.simulation_downlink.power_control()
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[0], 
                            [40], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[1], 
                            [40], atol=1e-2)
        
        self.simulation_downlink.calculate_sinr()      
        npt.assert_allclose(self.simulation_downlink.receiver.rx_power,
                            [-70.45, -79.33], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.rx_interference,
                            [-90.27, -92.05], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.thermal_noise,
                            -111.96*np.ones(2), atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.total_interference,
                            [-90.24, -92.01], atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.sinr,
                            [19.79, 12.68], atol=1e-2)
        
    def test_simulation_2bs_4ue(self):
        self.param.num_base_stations = 1
        self.param.num_clusters = 2
        self.param.ue_k = 2
        self.param.ue_k_m = 1
        self.param.mcl = 122
        self.param.ho_margin = 3
        self.simulation_downlink = SimulationDownlink(self.param)
        
        # after object instatiation, transmitter and receiver are only arrays
        self.assertEqual(len(self.simulation_downlink.transmitter), 2)
        self.assertEqual(len(self.simulation_downlink.receiver), 4)
        self.assertEqual(self.simulation_downlink.coupling_loss.shape, (2,4))
        
        # after initialize(), transmitters must be created
        self.simulation_downlink.initialize()
        self.assertEqual(self.simulation_downlink.transmitter.num_stations, 2)
        
        # it is time to create user equipments
        self.simulation_downlink.create_ue()
        self.simulation_downlink.receiver.x = np.array([-2000, -500, 400, 1500])
        self.simulation_downlink.receiver.y = np.array([0, 0, 0, 0])
        self.assertEqual(self.simulation_downlink.receiver.num_stations, 4)

        self.simulation_downlink.transmitter.tx_antenna = [Antenna(0), Antenna(1)]
        self.simulation_downlink.receiver.rx_antenna = [Antenna(2), Antenna(3), Antenna(4), Antenna(5)]
        
        # let's calculate coupling loss
        self.simulation_downlink.calculate_coupling_loss()
        npt.assert_allclose(self.simulation_downlink.coupling_loss, 
                            np.array([[124.42,  117.40,  125.35,  129.38], [132.97,  125.95,  116.99,  114.40]]), 
                            atol=1e-2)
        
        # Now we connect base stations to user equipments
        self.simulation_downlink.connect_ue_to_bs()
        self.assertEqual(self.simulation_downlink.link, {0: [0,1], 1: [2,3]})
        
        # Scheduling algorirhm
        self.simulation_downlink.scheduler()
        npt.assert_equal(self.simulation_downlink.receiver.bandwidth,
                         100*np.ones(4))
        
        # apply power control to set transmit powers
        self.simulation_downlink.power_control()
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[0], 
                            [36.98,36.98], 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.transmitter.tx_power[1], 
                            [36.98,36.98], 
                            atol=1e-2)        
        
        # calculate SINR
        self.simulation_downlink.calculate_sinr()
        npt.assert_allclose(self.simulation_downlink.receiver.rx_power,
                            [-87.43, -80.41, -80.00, -77.41],
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.rx_interference,
                            [-92.97,  -85.95,  -85.35, -89.38], 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.thermal_noise,
                            -114.97*np.ones(4), 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.total_interference,
                            [-92.94, -85.94, -85.34, -89.37], 
                            atol=1e-2)
        npt.assert_allclose(self.simulation_downlink.receiver.sinr,
                            [5.50, 5.52, 5.34, 11.95], 
                            atol=1e-2)
        
if __name__ == '__main__':
    unittest.main()
                