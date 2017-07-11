# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 18:32:30 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt
import math

from sharc.simulation_uplink import SimulationUplink
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.station_factory import StationFactory

class SimulationUplinkTest(unittest.TestCase):
    
    def setUp(self):
        self.param = ParametersImt()
        self.param.topology = "SINGLE_BS"
        self.param.num_macrocell_sites = 19
        self.param.num_clusters = 2
        self.param.intersite_distance = 200
        self.param.minimum_separation_distance_bs_ue = 10
        self.param.interfered_with = False
        self.param.frequency = 10000
        self.param.bandwidth = 100
        self.param.rb_bandwidth = 0.180
        self.param.guard_band_ratio = 0.1
        self.param.ho_margin = 3
        self.param.bs_load_probability = 1
        self.param.num_resource_blocks = 10
        self.param.bs_tx_power = 40
        self.param.bs_height = 6
        self.param.bs_aclr = 40
        self.param.bs_acs = 30
        self.param.bs_noise_figure = 7
        self.param.bs_noise_temperature = 290
        self.param.bs_feed_loss = 3
        self.param.ul_attenuation_factor = 0.4
        self.param.ul_sinr_min = -10
        self.param.ul_sinr_max = 22
        self.param.ue_k = 2
        self.param.ue_k_m = 1
        self.param.ue_tx_power_control = "OFF"
        self.param.ue_tx_power_target = -95
        self.param.ue_tx_power_alfa = 0.8
        self.param.ue_tx_power = 20
        self.param.ue_height = 1.5
        self.param.ue_aclr = 35
        self.param.ue_acs = 25
        self.param.ue_noise_figure = 9
        self.param.ue_feed_loss = 3
        self.param.dl_attenuation_factor = 0.6
        self.param.dl_sinr_min = -10
        self.param.dl_sinr_max = 30
        self.param.channel_model = "FSPL"
        self.param.line_of_sight_prob = 0.75 # probability of line-of-sight (not for FSPL)
        self.param.noise_temperature = 290
        self.param.BOLTZMANN_CONSTANT = 1.38064852e-23
        
        self.param_ant = ParametersAntennaImt()
        self.param_ant.bs_tx_element_max_g = 5
        self.param_ant.bs_tx_element_phi_3db = 80
        self.param_ant.bs_tx_element_theta_3db = 65
        self.param_ant.bs_tx_element_am = 30
        self.param_ant.bs_tx_element_sla_v = 30
        self.param_ant.bs_tx_n_rows = 8
        self.param_ant.bs_tx_n_columns = 8
        self.param_ant.bs_tx_element_horiz_spacing = 0.5
        self.param_ant.bs_tx_element_vert_spacing = 0.5
        self.param_ant.bs_rx_element_max_g = 5
        self.param_ant.bs_rx_element_phi_3db = 65
        self.param_ant.bs_rx_element_theta_3db = 65
        self.param_ant.bs_rx_element_am = 30
        self.param_ant.bs_rx_element_sla_v = 30
        self.param_ant.bs_rx_n_rows = 8
        self.param_ant.bs_rx_n_columns = 8
        self.param_ant.bs_rx_element_horiz_spacing = 0.5
        self.param_ant.bs_rx_element_vert_spacing = 0.5
        self.param_ant.ue_tx_element_max_g = 5
        self.param_ant.ue_tx_element_phi_3db = 80
        self.param_ant.ue_tx_element_theta_3db = 65
        self.param_ant.ue_tx_element_am = 30
        self.param_ant.ue_tx_element_sla_v = 30
        self.param_ant.ue_tx_n_rows = 4
        self.param_ant.ue_tx_n_columns = 4
        self.param_ant.ue_tx_element_horiz_spacing = 0.5
        self.param_ant.ue_tx_element_vert_spacing = 0.5
        self.param_ant.ue_rx_element_max_g = 5
        self.param_ant.ue_rx_element_phi_3db = 90
        self.param_ant.ue_rx_element_theta_3db = 90
        self.param_ant.ue_rx_element_am = 25
        self.param_ant.ue_rx_element_sla_v = 25
        self.param_ant.ue_rx_n_rows = 4
        self.param_ant.ue_rx_n_columns = 4
        self.param_ant.ue_rx_element_horiz_spacing = 0.5
        self.param_ant.ue_rx_element_vert_spacing = 0.5
        
        self.simulation = SimulationUplink(self.param, ParametersFss(), self.param_ant)
        self.simulation.initialize()
        
        
    def test_simulation_2bs_4ue(self):
        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param,
                                                                       self.param_ant,
                                                                       self.simulation.topology)
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)
        
        self.simulation.ue = StationFactory.generate_imt_ue(self.param,
                                                            self.param_ant,
                                                            self.simulation.topology)
        self.simulation.ue.x = np.array([20, 70, 110, 170])
        self.simulation.ue.y = np.array([ 0,  0,   0,   0])
        self.simulation.ue.antenna = np.array([AntennaOmni(10), AntennaOmni(11), AntennaOmni(22), AntennaOmni(23)])
        self.simulation.ue.active = np.ones(4, dtype=bool)
        
        # test connection method
        self.simulation.connect_ue_to_bs()
        self.assertEqual(self.simulation.link, {0: [0,1], 1: [2,3]})
        
        # We do not test the selection method here because in this specific 
        # scenario we do not want to change the order of the UE's 
        #self.simulation.select_ue()
        
        # test coupling loss method
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs, 
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        npt.assert_allclose(self.simulation.coupling_loss_imt, 
                            np.array([[78.47-1-10,  89.35-1-11,  93.27-1-22,  97.05-1-23], 
                                      [97.55-2-10,  94.72-2-11,  91.53-2-22,  81.99-2-23]]), 
                            atol=1e-2)
        
        # test scheduler and bandwidth allocation
        self.simulation.scheduler()
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)       
        npt.assert_allclose(self.simulation.ue.bandwidth, bandwidth_per_ue*np.ones(4), atol=1e-2)
        
        # there is no power control, so UE's will transmit at maximum power
        self.simulation.power_control()
        npt.assert_allclose(self.simulation.ue.tx_power, 20*np.ones(4))
        
        # test method that calculates SINR 
        self.simulation.calculate_sinr()
        # check BS received power
        npt.assert_allclose(self.simulation.bs.rx_power[0], 
                            np.array([20-(78.47-1-10), 20-(89.35-1-11)]),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_power[1], 
                            np.array([20-(91.53-2-22), 20-(81.99-2-23)]),
                            atol=1e-2)        
        # check BS received interference
        npt.assert_allclose(self.simulation.bs.rx_interference[0], 
                            np.array([20-(93.27-1-22),  20-(97.05-1-23)]),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_interference[1], 
                            np.array([20-(97.55-2-10), 20-(94.72-2-11)]),
                            atol=1e-2)      
        # check BS thermal noise
        npt.assert_allclose(self.simulation.bs.thermal_noise, 
                            10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e6) + 7,
                            atol=1e-2)
        # check BS thermal noise + interference
        npt.assert_allclose(self.simulation.bs.total_interference[0], 
                            10*np.log10(np.power(10, 0.1*np.array([20-(93.27-1-22),  20-(97.05-1-23)])) +
                                        np.power(10, 0.1*(-120.44))),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.total_interference[1], 
                            10*np.log10(np.power(10, 0.1*np.array([20-(97.55-2-10), 20-(94.72-2-11)])) +
                                        np.power(10, 0.1*(-120.44))),
                            atol=1e-2)    
        # check SNR 
        npt.assert_allclose(self.simulation.bs.snr[0], 
                            np.array([20-(78.47-1-10) - (-120.44),  20-(89.35-1-11) - (-120.44)]),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.snr[1], 
                            np.array([20-(91.53-2-22) - (-120.44),  20-(81.99-2-23) - (-120.44)]),
                            atol=1e-2)
        # check SINR (using only the I term since I >> N in this case)
        npt.assert_allclose(self.simulation.bs.sinr[0], 
                            np.array([20-(78.47-1-10) - (20-(93.27-1-22)),  20-(89.35-1-11) - (20-(97.05-1-23))]),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.sinr[1], 
                            np.array([20-(91.53-2-22) - (20-(97.55-2-10)),  20-(81.99-2-23) - (20-(94.72-2-11))]),
                            atol=1e-2)


    def test_calculate_gains(self):
        pass
    
    
    def test_calculate_coupling_loss(self):
        pass
    
    
    def test_calculate_imt_ul_tput(self):
        pass
    
    
    
                
if __name__ == '__main__':
    unittest.main()
    
