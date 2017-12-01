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
from sharc.parameters.parameters import Parameters
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.station_factory import StationFactory

class SimulationDownlinkHapsTest(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()
        
        self.param.general.imt_link = "DOWNLINK"
        
        self.param.imt.topology = "SINGLE_BS"
        self.param.imt.num_macrocell_sites = 19
        self.param.imt.num_clusters = 2
        self.param.imt.intersite_distance = 150
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 10000
        self.param.imt.bandwidth = 100
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.ho_margin = 3
        self.param.imt.bs_load_probability = 1
        self.param.imt.num_resource_blocks = 10
        self.param.imt.bs_conducted_power = 10
        self.param.imt.bs_height = 6
        self.param.imt.bs_aclr = 40
        self.param.imt.bs_acs = 30
        self.param.imt.bs_noise_figure = 7
        self.param.imt.bs_noise_temperature = 290
        self.param.imt.bs_feed_loss = 3
        self.param.imt.ul_attenuation_factor = 0.4
        self.param.imt.ul_sinr_min = -10
        self.param.imt.ul_sinr_max = 22
        self.param.imt.ue_k = 2
        self.param.imt.ue_k_m = 1
        self.param.imt.ue_indoor_percent = 0
        self.param.imt.ue_distribution_distance = "RAYLEIGH"
        self.param.imt.ue_distribution_azimuth = "UNIFORM"
        self.param.imt.ue_tx_power_control = "OFF"
        self.param.imt.ue_p_o_pusch = -95
        self.param.imt.ue_alfa = 0.8
        self.param.imt.ue_p_cmax = 20
        self.param.imt.ue_conducted_power = 10
        self.param.imt.ue_height = 1.5
        self.param.imt.ue_aclr = 35
        self.param.imt.ue_acs = 25
        self.param.imt.ue_noise_figure = 9
        self.param.imt.ue_feed_loss = 3
        self.param.imt.ue_body_loss = 4
        self.param.imt.dl_attenuation_factor = 0.6
        self.param.imt.dl_sinr_min = -10
        self.param.imt.dl_sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        self.param.imt.line_of_sight_prob = 0.75 # probability of line-of-sight (not for FSPL)
        self.param.imt.shadowing = False
        self.param.imt.noise_temperature = 290
        self.param.imt.BOLTZMANN_CONSTANT = 1.38064852e-23
        
        self.param.antenna_imt.bs_tx_element_max_g = 10
        self.param.antenna_imt.bs_tx_element_phi_3db = 80
        self.param.antenna_imt.bs_tx_element_theta_3db = 80
        self.param.antenna_imt.bs_tx_element_am = 25
        self.param.antenna_imt.bs_tx_element_sla_v = 25
        self.param.antenna_imt.bs_tx_n_rows = 16
        self.param.antenna_imt.bs_tx_n_columns = 16
        self.param.antenna_imt.bs_tx_element_horiz_spacing = 1
        self.param.antenna_imt.bs_tx_element_vert_spacing = 1
        self.param.antenna_imt.bs_rx_element_max_g = 5
        self.param.antenna_imt.bs_rx_element_phi_3db = 65
        self.param.antenna_imt.bs_rx_element_theta_3db = 65
        self.param.antenna_imt.bs_rx_element_am = 30
        self.param.antenna_imt.bs_rx_element_sla_v = 30
        self.param.antenna_imt.bs_rx_n_rows = 2
        self.param.antenna_imt.bs_rx_n_columns = 2
        self.param.antenna_imt.bs_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_rx_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_tx_element_max_g = 5
        self.param.antenna_imt.ue_tx_element_phi_3db = 65
        self.param.antenna_imt.ue_tx_element_theta_3db = 65
        self.param.antenna_imt.ue_tx_element_am = 30
        self.param.antenna_imt.ue_tx_element_sla_v = 30
        self.param.antenna_imt.ue_tx_n_rows = 2
        self.param.antenna_imt.ue_tx_n_columns = 1
        self.param.antenna_imt.ue_tx_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_tx_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_rx_element_max_g = 10
        self.param.antenna_imt.ue_rx_element_phi_3db = 90
        self.param.antenna_imt.ue_rx_element_theta_3db = 90
        self.param.antenna_imt.ue_rx_element_am = 25
        self.param.antenna_imt.ue_rx_element_sla_v = 25
        self.param.antenna_imt.ue_rx_n_rows = 16
        self.param.antenna_imt.ue_rx_n_columns = 16
        self.param.antenna_imt.ue_rx_element_horiz_spacing = 1
        self.param.antenna_imt.ue_rx_element_vert_spacing = 1
        
        self.param.haps.frequency = 10000
        self.param.haps.bandwidth = 200
        self.param.haps.altitude = 20000
        self.param.haps.lat_deg = 0
        self.param.haps.elevation = 270
        self.param.haps.azimuth = 0
        self.param.haps.eirp_density = 4.4
        self.param.haps.inr_scaling = 1
        self.param.haps.antenna_gain = 28
        self.param.haps.tx_power_density = self.param.haps.eirp_density - self.param.haps.antenna_gain - 60
        self.param.haps.antenna_pattern = "OMNI"
        self.param.haps.imt_altitude = 0
        self.param.haps.imt_lat_deg = 0
        self.param.haps.imt_long_diff_deg = 0
        self.param.haps.season = "SUMMER"
        self.param.haps.channel_model = "FSPL"
        self.param.haps.antenna_l_n = -25
        self.param.haps.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.haps.EARTH_RADIUS = 6371000
        
        

    def test_simulation_1bs_1haps(self):
        """
        Test the interference generated by one HAPS (airbone) station to
        one IMT base station
        """
        self.param.general.system = "HAPS"
        
        self.simulation = SimulationDownlink(self.param)
        self.simulation.initialize()
        
        
        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0
        
        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology)
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)
        
        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
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
        
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2) 
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9
        
        self.simulation.system = StationFactory.generate_haps(self.param.haps, 0)
        
        # now we evaluate interference from HAPS to IMT UE
        self.simulation.calculate_sinr_ext()
        npt.assert_allclose(self.simulation.coupling_loss_imt_system, 
                            np.array([138.47-28-10,  138.47-28-11,  138.47-28-22,  138.47-28-23]), 
                            atol=1e-2)

        system_tx_power = (4.4 - 28 - 60) + 10*math.log10(bandwidth_per_ue*1e6) + 30
        npt.assert_allclose(self.simulation.ue.ext_interference, 
                            np.array([system_tx_power - (138.47-28-10) - 7,  system_tx_power - (138.47-28-11) - 7,  system_tx_power - (138.47-28-22) - 7,  system_tx_power - (138.47-28-23) - 7]), 
                            atol=1e-2)
        
        interference = 10*np.log10(np.power(10, 0.1*np.array([ -86.71, -85.06, -76.03, -78.59 ])) \
                                 + np.power(10, 0.1*np.array([ -84.53, -83.53, -72.53, -71.53 ])))
        
        npt.assert_allclose(self.simulation.ue.sinr_ext, 
                            np.array([-73.48, -83.36, -73.54, -63.00]) - interference, 
                            atol=1e-2)       
        
        npt.assert_allclose(self.simulation.ue.inr, 
                            np.array([-84.53, -83.53, -72.53, -71.53]) - thermal_noise, 
                            atol=1e-2)
    
        
        
        
        
if __name__ == '__main__':
    unittest.main()
