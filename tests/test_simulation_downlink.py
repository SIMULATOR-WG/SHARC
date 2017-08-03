# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 16:22:42 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt
import math

from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.station_factory import StationFactory

class SimulationDownlinkTest(unittest.TestCase):

    def setUp(self):
        self.param = ParametersImt()
        self.param.topology = "SINGLE_BS"
        self.param.num_macrocell_sites = 19
        self.param.num_clusters = 2
        self.param.intersite_distance = 150
        self.param.minimum_separation_distance_bs_ue = 10
        self.param.interfered_with = False
        self.param.frequency = 10000
        self.param.bandwidth = 100
        self.param.rb_bandwidth = 0.180
        self.param.guard_band_ratio = 0.1
        self.param.ho_margin = 3
        self.param.bs_load_probability = 1
        self.param.num_resource_blocks = 10
        self.param.bs_conducted_power = 10
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
        self.param.ue_indoor_percent = 0
        self.param.ue_distribution_distance = "RAYLEIGH"
        self.param.ue_distribution_azimuth = "UNIFORM"
        self.param.ue_tx_power_control = "OFF"
        self.param.ue_p_o_pusch = -95
        self.param.ue_alfa = 0.8
        self.param.ue_p_cmax = 20
        self.param.ue_conducted_power = 10
        self.param.ue_height = 1.5
        self.param.ue_aclr = 35
        self.param.ue_acs = 25
        self.param.ue_noise_figure = 9
        self.param.ue_feed_loss = 3
        self.param.ue_body_loss = 4
        self.param.dl_attenuation_factor = 0.6
        self.param.dl_sinr_min = -10
        self.param.dl_sinr_max = 30
        self.param.channel_model = "FSPL"
        self.param.line_of_sight_prob = 0.75 # probability of line-of-sight (not for FSPL)
        self.param.noise_temperature = 290
        self.param.BOLTZMANN_CONSTANT = 1.38064852e-23
        
        self.param_ant = ParametersAntennaImt()
        self.param_ant.bs_tx_element_max_g = 10
        self.param_ant.bs_tx_element_phi_3db = 80
        self.param_ant.bs_tx_element_theta_3db = 80
        self.param_ant.bs_tx_element_am = 25
        self.param_ant.bs_tx_element_sla_v = 25
        self.param_ant.bs_tx_n_rows = 16
        self.param_ant.bs_tx_n_columns = 16
        self.param_ant.bs_tx_element_horiz_spacing = 1
        self.param_ant.bs_tx_element_vert_spacing = 1
        self.param_ant.bs_rx_element_max_g = 5
        self.param_ant.bs_rx_element_phi_3db = 65
        self.param_ant.bs_rx_element_theta_3db = 65
        self.param_ant.bs_rx_element_am = 30
        self.param_ant.bs_rx_element_sla_v = 30
        self.param_ant.bs_rx_n_rows = 2
        self.param_ant.bs_rx_n_columns = 2
        self.param_ant.bs_rx_element_horiz_spacing = 0.5
        self.param_ant.bs_rx_element_vert_spacing = 0.5
        self.param_ant.ue_tx_element_max_g = 5
        self.param_ant.ue_tx_element_phi_3db = 65
        self.param_ant.ue_tx_element_theta_3db = 65
        self.param_ant.ue_tx_element_am = 30
        self.param_ant.ue_tx_element_sla_v = 30
        self.param_ant.ue_tx_n_rows = 2
        self.param_ant.ue_tx_n_columns = 1
        self.param_ant.ue_tx_element_horiz_spacing = 0.5
        self.param_ant.ue_tx_element_vert_spacing = 0.5
        self.param_ant.ue_rx_element_max_g = 10
        self.param_ant.ue_rx_element_phi_3db = 90
        self.param_ant.ue_rx_element_theta_3db = 90
        self.param_ant.ue_rx_element_am = 25
        self.param_ant.ue_rx_element_sla_v = 25
        self.param_ant.ue_rx_n_rows = 16
        self.param_ant.ue_rx_n_columns = 16
        self.param_ant.ue_rx_element_horiz_spacing = 1
        self.param_ant.ue_rx_element_vert_spacing = 1
        
        self.param_fss_ss = ParametersFss()
        self.param_fss_ss.frequency = 10000
        self.param_fss_ss.bandwidth = 100
        self.param_fss_ss.sat_altitude = 35786000
        self.param_fss_ss.sat_lat_deg = 0
        self.param_fss_ss.sat_noise_temperature = 950
        self.param_fss_ss.sat_interference_noise_ratio = -12.2
        self.param_fss_ss.sat_rx_antenna_gain = 51
        self.param_fss_ss.sat_rx_antenna_pattern = "OMNI"
        self.param_fss_ss.imt_altitude = 1000
        self.param_fss_ss.imt_lat_deg = -23.5629739
        self.param_fss_ss.imt_long_diff_deg = (-46.6555132-75)
        self.param_fss_ss.channel_model = "FSPL"
        self.param_fss_ss.line_of_sight_prob = 0.01
        self.param_fss_ss.surf_water_vapour_density = 7.5
        self.param_fss_ss.specific_gaseous_att = 0.1
        self.param_fss_ss.time_ratio = 0.5
        self.param_fss_ss.sat_rx_antenna_l_s = -20    
        self.param_fss_ss.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param_fss_ss.EARTH_RADIUS = 6371000        
        
        self.param_fss_es = ParametersFssEs()
        self.param_fss_es.x = -5000
        self.param_fss_es.y = 0
        self.param_fss_es.height = 10
        self.param_fss_es.elevation = 20
        self.param_fss_es.azimuth = 0
        self.param_fss_es.frequency = 10000
        self.param_fss_es.bandwidth = 100
        self.param_fss_es.tx_power_density = -60
        self.param_fss_es.antenna_gain = 50
        self.param_fss_es.antenna_pattern = "OMNI"
        self.param_fss_es.channel_model = "FSPL"
        self.param_fss_es.line_of_sight_prob = 1 
        self.param_fss_es.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param_fss_es.EARTH_RADIUS = 6371000        


    def test_simulation_2bs_4ue_fss_ss(self):
        self.simulation = SimulationDownlink(self.param, self.param_fss_ss, self.param_ant)
        self.simulation.initialize()

        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0
        
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
        
        # there is no power control, so BS's will transmit at maximum power
        self.simulation.power_control()
        p_tx = 10 + 0 - 3 - 10*math.log10(2)
        npt.assert_allclose(self.simulation.bs.tx_power[0], np.array([p_tx, p_tx]), atol=1e-2)
        npt.assert_allclose(self.simulation.bs.tx_power[1], np.array([p_tx, p_tx]), atol=1e-2)
        
        # test method that calculates SINR 
        self.simulation.calculate_sinr()
        # check UE received power
        npt.assert_allclose(self.simulation.ue.rx_power, 
                            np.array([p_tx-(78.47-1-10)-7, p_tx-(89.35-1-11)-7, p_tx-(91.53-2-22)-7, p_tx-(81.99-2-23)-7]),
                            atol=1e-2)
        # check UE received interference
        npt.assert_allclose(self.simulation.ue.rx_interference, 
                            np.array([p_tx-(97.55-2-10)-7,  p_tx-(94.72-2-11)-7, p_tx-(93.27-1-22)-7, p_tx-(97.05-1-23)-7]),
                            atol=1e-2)
        # check UE thermal noise
        npt.assert_allclose(self.simulation.ue.thermal_noise, 
                            10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9,
                            atol=1e-2)
        # check BS thermal noise + interference
        npt.assert_allclose(self.simulation.ue.total_interference, 
                            10*np.log10(np.power(10, 0.1*np.array([p_tx-(97.55-2-10)-7,  p_tx-(94.72-2-11)-7, p_tx-(93.27-1-22)-7, p_tx-(97.05-1-23)-7])) +
                                        np.power(10, 0.1*(-88.44))),
                            atol=1e-2)
        # check SNR 
        npt.assert_allclose(self.simulation.ue.snr, 
                            np.array([-70.48 - (-88.44),  -80.36 - (-88.44), -70.54 - (-88.44),  -60.00 - (-88.44)]),
                            atol=1e-2)        
        # check SINR
        npt.assert_allclose(self.simulation.ue.sinr, 
                            np.array([-70.48 - (-85.49), -80.36 - (-83.19), -70.54 - (-73.15), -60.00 - (-75.82)]),
                            atol=1e-2)        

        self.simulation.system = StationFactory.generate_fss_space_station(self.param_fss_ss)
        self.simulation.system.x = np.array([0])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param_fss_ss.sat_altitude])
        
        # test the method that calculates interference from IMT UE to FSS space station
        self.simulation.calculate_external_interference()
        # check coupling loss
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            np.array([203.52-51-1, 203.52-51-1, 203.52-51-2, 203.52-51-2]),
                            atol=1e-2)
        # check interference generated by BS to FSS space station
        interference = 10 - 10*np.log10(2) - np.array([203.52-51-1, 203.52-51-1, 203.52-51-2, 203.52-51-2])- 3 + 10*math.log10(45/100)
        rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference)))
        self.assertAlmostEqual(self.simulation.system.rx_interference,
                               rx_interference,
                               delta=.01)
        # check FSS space station thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*950*1e3*100*1e6)
        self.assertAlmostEqual(self.simulation.system.thermal_noise, 
                               thermal_noise,
                               delta=.01)      
        # check INR at FSS space station
        self.assertAlmostEqual(self.simulation.system.inr, 
                               -144.448 - (-88.821),
                               delta=.01)        
       
        
    def test_simulation_2bs_4ue_fss_es(self):
        self.simulation = SimulationDownlink(self.param, self.param_fss_es, self.param_ant)
        self.simulation.initialize()
        
        
        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0
        
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
        
        self.simulation.connect_ue_to_bs()
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs, 
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        self.simulation.scheduler()
        self.simulation.power_control()
        self.simulation.calculate_sinr()
        # check SINR
        npt.assert_allclose(self.simulation.ue.sinr, 
                            np.array([-70.48 - (-85.49), -80.36 - (-83.19), -70.54 - (-73.15), -60.00 - (-75.82)]),
                            atol=1e-2)        

        self.simulation.system = StationFactory.generate_fss_earth_station(self.param_fss_es)
        self.simulation.system.x = np.array([-2000])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param_fss_es.height])
        
        self.simulation.calculate_sinr_ext()
        npt.assert_allclose(self.simulation.coupling_loss_imt_system, 
                            np.array([118.55-50-10,  118.76-50-11,  118.93-50-22,  119.17-50-23]), 
                            atol=1e-2)

        bw = 100*1e6*0.9/2
        system_tx_power = -60 + 10*math.log10(bw) + 30
        npt.assert_allclose(self.simulation.ue.ext_interference, 
                            np.array([system_tx_power - (118.55-50-10) - 7,  system_tx_power - (118.76-50-11) - 7,  system_tx_power - (118.93-50-22) - 7,  system_tx_power - (119.17-50-23) - 7]), 
                            atol=1e-2)
        
        interference = 10*np.log10(np.power(10, 0.1*np.array([ -85.49, -83.19, -73.15, -75.82 ])) \
                                 + np.power(10, 0.1*np.array([ -19.02, -18.23,  -7.40,  -6.64 ])))
        
        npt.assert_allclose(self.simulation.ue.sinr_ext, 
                            np.array([-70.48, -80.36, -70.54, -60.00]) - interference, 
                            atol=1e-2)       
        
        
if __name__ == '__main__':
    unittest.main()
