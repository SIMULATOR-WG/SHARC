# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 15:43:06 2018

@author: Calil
"""


import unittest
import numpy as np
import numpy.testing as npt
import os.path as path

from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters import Parameters
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.station_factory import StationFactory
from sharc.propagation.propagation_factory import PropagationFactory

class SimulationIndoorTest(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()

        self.param.general.imt_link = "DOWNLINK"
        self.param.general.enable_cochannel = True
        self.param.general.enable_adjacent_channel = False
        self.param.general.overwrite_output = True

        self.param.imt.topology = "INDOOR"
        self.param.imt.num_macrocell_sites = 19
        self.param.imt.num_clusters = 1
        self.param.imt.intersite_distance = 339
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 40000
        self.param.imt.bandwidth = 200
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "ITU 265-E"
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.bs_load_probability = 1
        self.param.imt.num_resource_blocks = 10
        self.param.imt.bs_conducted_power = 2
        self.param.imt.bs_height = 3
        self.param.imt.bs_noise_figure = 12
        self.param.imt.bs_noise_temperature = 290
        self.param.imt.bs_ohmic_loss = 3
        self.param.imt.ul_attenuation_factor = 0.4
        self.param.imt.ul_sinr_min = -10
        self.param.imt.ul_sinr_max = 22
        self.param.imt.ue_k = 1
        self.param.imt.ue_k_m = 1
        self.param.imt.ue_indoor_percent = 95
        self.param.imt.ue_distribution_type = "ANGLE_AND_DISTANCE"
        self.param.imt.ue_distribution_distance = "RAYLEIGH"
        self.param.imt.ue_distribution_azimuth = "UNIFORM"
        self.param.imt.ue_tx_power_control = "OFF"
        self.param.imt.ue_p_o_pusch = -95
        self.param.imt.ue_alpha = 1
        self.param.imt.ue_p_cmax = 22
        self.param.imt.ue_height = 1.5
        self.param.imt.ue_noise_figure = 12
        self.param.imt.ue_ohmic_loss = 3
        self.param.imt.ue_body_loss = 4
        self.param.imt.dl_attenuation_factor = 0.6
        self.param.imt.dl_sinr_min = -10
        self.param.imt.dl_sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        self.param.imt.shadowing = False
        self.param.imt.noise_temperature = 290
        self.param.imt.BOLTZMANN_CONSTANT = 1.38064852e-23

        self.param.antenna_imt.normalization = False
        self.param.antenna_imt.bs_normalization_file = path.join('..','sharc','antenna','beamforming_normalization','bs_indoor_norm.npz')
        self.param.antenna_imt.ue_normalization_file = path.join('..','sharc','antenna','beamforming_normalization','ue_norm.npz')
        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_tx_element_max_g = 5
        self.param.antenna_imt.bs_tx_element_phi_3db = 90
        self.param.antenna_imt.bs_tx_element_theta_3db = 90
        self.param.antenna_imt.bs_tx_element_am = 25
        self.param.antenna_imt.bs_tx_element_sla_v = 25
        self.param.antenna_imt.bs_tx_n_rows = 8
        self.param.antenna_imt.bs_tx_n_columns = 16
        self.param.antenna_imt.bs_tx_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_tx_element_vert_spacing = 0.5
        self.param.antenna_imt.bs_rx_element_max_g = 5
        self.param.antenna_imt.bs_rx_element_phi_deg_3db = 90
        self.param.antenna_imt.bs_rx_element_theta_deg_3db = 90
        self.param.antenna_imt.bs_rx_element_am = 25
        self.param.antenna_imt.bs_rx_element_sla_v = 25
        self.param.antenna_imt.bs_rx_n_rows = 8
        self.param.antenna_imt.bs_rx_n_columns = 16
        self.param.antenna_imt.bs_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_rx_element_vert_spacing = 0.5
        self.param.antenna_imt.bs_downtilt_deg = 10
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_tx_element_max_g = 5
        self.param.antenna_imt.ue_tx_element_phi_deg_3db = 90
        self.param.antenna_imt.ue_tx_element_theta_deg_3db = 90
        self.param.antenna_imt.ue_tx_element_am = 25
        self.param.antenna_imt.ue_tx_element_sla_v = 25
        self.param.antenna_imt.ue_tx_n_rows = 4
        self.param.antenna_imt.ue_tx_n_columns = 4
        self.param.antenna_imt.ue_tx_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_tx_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_rx_element_max_g = 5
        self.param.antenna_imt.ue_rx_element_phi_3db = 90
        self.param.antenna_imt.ue_rx_element_theta_3db = 90
        self.param.antenna_imt.ue_rx_element_am = 25
        self.param.antenna_imt.ue_rx_element_sla_v = 25
        self.param.antenna_imt.ue_rx_n_rows = 4
        self.param.antenna_imt.ue_rx_n_columns = 4
        self.param.antenna_imt.ue_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_rx_element_vert_spacing = 0.5

        self.param.indoor.basic_path_loss = "FSPL"
        self.param.indoor.n_rows = 1
        self.param.indoor.n_colums = 1
        self.param.indoor.street_width = 30
        self.param.indoor.ue_indoor_percent = 0.95
        self.param.indoor.building_class = "TRADITIONAL"
        self.param.indoor.intersite_distance = 30
        self.param.indoor.num_cells = 4
        self.param.indoor.num_floors = 1

        self.param.fss_es.x = 135
        self.param.fss_es.y = 65
        self.param.fss_es.location = "FIXED"
        self.param.fss_es.height = 10
        self.param.fss_es.elevation_min = 10
        self.param.fss_es.elevation_max = 10
        self.param.fss_es.azimuth = "-180"
        self.param.fss_es.frequency = 40000
        self.param.fss_es.bandwidth = 180
        self.param.fss_es.noise_temperature = 400
        self.param.fss_es.tx_power_density = -69
        self.param.fss_es.antenna_gain = 47
        self.param.fss_es.antenna_pattern = "ITU-R S.580"
        self.param.fss_es.channel_model = "FSPL"
        self.param.fss_es.line_of_sight_prob = 1
        self.param.fss_es.adjacent_ch_selectivity = 0
        self.param.fss_es.diameter = 0.74
        self.param.fss_es.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.fss_es.EARTH_RADIUS = 6371000


    def test_simulation_fss_es(self):
        # Initialize stations
        self.param.general.system = "FSS_ES"

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState(101)

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)
        self.assertTrue(np.all(self.simulation.bs.active))
        
        self.simulation.system = StationFactory.generate_fss_earth_station(self.param.fss_es,
                                                                           random_number_gen)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        
#        print("Random position:")
#        self.simulation.plot_scenario()
        self.simulation.ue.x = np.array([0.0, 45.0, 75.0,120.0])
        self.simulation.ue.y = np.array([0.0, 50.0,  0.0, 50.0])
#        print("Forced position:")
#        self.simulation.plot_scenario()
        
        # Connect and select UEs
        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        self.assertTrue(np.all(self.simulation.ue.active))
        self.assertDictEqual(self.simulation.link,{0:[0],1:[1],2:[2],3:[3]})
        
        # Test BS-to-UE angles in the IMT coord system
        expected_azi = np.array([[-120.96,  39.80, -22.62,  13.39],
                                 [-150.95,  90.00, -39.81,  18.43],
                                 [-161.57, 140.19, -90.00,  29.06],
                                 [-166.61, 157.38,-140.19,  59.03]])
        npt.assert_allclose(self.simulation.bs_to_ue_phi,
                            expected_azi,
                            atol=1e-2)
        expected_ele = np.array([[92.95,  92.20, 91.32,  90.79],
                                 [91.67,  93.43, 92.20,  91.09],
                                 [91.09,  92.20, 93.43,  91.67],
                                 [90.79,  91.32, 92.20,  92.95]])
        npt.assert_allclose(self.simulation.bs_to_ue_theta,
                            expected_ele,
                            atol=1e-2)
        
        # Test BS-to-UE angles in the local coord system
        expected_loc = [(np.array([-86.57]),np.array([120.92])),
                        (np.array([ 86.57]),np.array([ 90.00])),
                        (np.array([-86.57]),np.array([ 90.00])),
                        (np.array([ 86.57]),np.array([ 59.08]))]
        expected_beam = [(-86.57,30.92),
                         ( 86.57, 0.00),
                         (-86.57, 0.00),
                         ( 86.57,-30.92)]
        for k in range(self.simulation.bs.num_stations):
            
            self.assertEqual(self.simulation.bs.antenna[k].azimuth,0.0)
            self.assertEqual(self.simulation.bs.antenna[k].elevation,-90.0)
            
            lo_angles = self.simulation.bs.antenna[k].to_local_coord(expected_azi[k,k],
                                                                  expected_ele[k,k])
            npt.assert_array_almost_equal(lo_angles,expected_loc[k],decimal=2)
            npt.assert_array_almost_equal(self.simulation.bs.antenna[k].beams_list[0],
                                          expected_beam[k],decimal=2)
            
        # Test angle to ES in the IMT coord system
        phi_es, theta_es = self.simulation.bs.get_pointing_vector_to(self.simulation.system)
        expected_phi_es = np.array([[18.44],[23.96],[33.69],[53.13]])
        npt.assert_array_almost_equal(phi_es,expected_phi_es,decimal=2)
        expected_theta_es = np.array([[86.83],[85.94],[84.46],[82.03]])
        npt.assert_array_almost_equal(theta_es,expected_theta_es,decimal=2)
        
        # Test angle to ES in the local coord system
        expected_es_loc = [(np.array([99.92]),np.array([18.70])),
                           (np.array([99.92]),np.array([24.28])),
                           (np.array([99.92]),np.array([34.09])),
                           (np.array([99.92]),np.array([53.54]))]
        for k in range(self.simulation.bs.num_stations):
            lo_angles = self.simulation.bs.antenna[k].to_local_coord(expected_phi_es[k],
                                                                     expected_theta_es[k])
            npt.assert_array_almost_equal(lo_angles,expected_es_loc[k],decimal=2)
        
        # Test gain to ES
        calc_gain = self.simulation.calculate_gains(self.simulation.bs,
                                               self.simulation.system)
        for k in range(self.simulation.bs.num_stations):
            beam = 0
            exp_gain = self.simulation.bs.antenna[k]._beam_gain(expected_es_loc[k][0],
                                                                expected_es_loc[k][1],
                                                                beam)
            self.assertAlmostEqual(np.asscalar(calc_gain[k]),np.asscalar(exp_gain),places=1)
            
if __name__ == '__main__':
    unittest.main()