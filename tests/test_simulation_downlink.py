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
from sharc.propagation.propagation_factory import PropagationFactory

class SimulationDownlinkTest(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()

        self.param.general.imt_link = "DOWNLINK"
        self.param.general.seed = 101
        self.param.general.enable_cochannel = True
        self.param.general.enable_adjacent_channel = False
        self.param.general.overwrite_output = True

        self.param.imt.topology = "SINGLE_BS"
        self.param.imt.wrap_around = False
        self.param.imt.num_clusters = 2
        self.param.imt.intersite_distance = 150
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 10000
        self.param.imt.bandwidth = 100
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "IMT-2020"
        self.param.imt.spurious_emissions = -13
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.ho_margin = 3
        self.param.imt.bs_load_probability = 1
        self.param.imt.num_resource_blocks = 10
        self.param.imt.bs_conducted_power = 10
        self.param.imt.bs_height = 6
        self.param.imt.bs_acs = 30
        self.param.imt.bs_noise_figure = 7
        self.param.imt.bs_noise_temperature = 290
        self.param.imt.bs_ohmic_loss = 3
        self.param.imt.ul_attenuation_factor = 0.4
        self.param.imt.ul_sinr_min = -10
        self.param.imt.ul_sinr_max = 22
        self.param.imt.ue_k = 2
        self.param.imt.ue_k_m = 1
        self.param.imt.ue_indoor_percent = 0
        self.param.imt.ue_distribution_distance = "RAYLEIGH"
        self.param.imt.ue_distribution_azimuth = "UNIFORM"
        self.param.imt.ue_distribution_type = "ANGLE_AND_DISTANCE"
        self.param.imt.ue_tx_power_control = "OFF"
        self.param.imt.ue_p_o_pusch = -95
        self.param.imt.ue_alpha = 0.8
        self.param.imt.ue_p_cmax = 20
        self.param.imt.ue_conducted_power = 10
        self.param.imt.ue_height = 1.5
        self.param.imt.ue_acs = 25
        self.param.imt.ue_noise_figure = 9
        self.param.imt.ue_ohmic_loss = 3
        self.param.imt.ue_body_loss = 4
        self.param.imt.dl_attenuation_factor = 0.6
        self.param.imt.dl_sinr_min = -10
        self.param.imt.dl_sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        self.param.imt.line_of_sight_prob = 0.75 # probability of line-of-sight (not for FSPL)
        self.param.imt.shadowing = False
        self.param.imt.noise_temperature = 290
        self.param.imt.BOLTZMANN_CONSTANT = 1.38064852e-23

        self.param.antenna_imt.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.antenna_imt.normalization = False
        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_normalization_file = None
        self.param.antenna_imt.bs_minimum_array_gain = -200
        self.param.antenna_imt.bs_element_max_g = 10
        self.param.antenna_imt.bs_element_phi_3db = 80
        self.param.antenna_imt.bs_element_theta_3db = 80
        self.param.antenna_imt.bs_element_am = 25
        self.param.antenna_imt.bs_element_sla_v = 25
        self.param.antenna_imt.bs_n_rows = 16
        self.param.antenna_imt.bs_n_columns = 16
        self.param.antenna_imt.bs_element_horiz_spacing = 1
        self.param.antenna_imt.bs_element_vert_spacing = 1
        self.param.antenna_imt.bs_multiplication_factor = 12
        self.param.antenna_imt.bs_downtilt = 10
        
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_normalization_file = None
        self.param.antenna_imt.ue_minimum_array_gain = -200
        self.param.antenna_imt.ue_element_max_g = 5
        self.param.antenna_imt.ue_element_phi_3db = 65
        self.param.antenna_imt.ue_element_theta_3db = 65
        self.param.antenna_imt.ue_element_am = 30
        self.param.antenna_imt.ue_element_sla_v = 30
        self.param.antenna_imt.ue_n_rows = 2
        self.param.antenna_imt.ue_n_columns = 1
        self.param.antenna_imt.ue_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_multiplication_factor = 12

        self.param.fss_ss.frequency = 10000
        self.param.fss_ss.bandwidth = 100
        self.param.fss_ss.acs = 0
        self.param.fss_ss.altitude = 35786000
        self.param.fss_ss.lat_deg = 0
        self.param.fss_ss.azimuth = 0
        self.param.fss_ss.elevation = 270
        self.param.fss_ss.tx_power_density = -30
        self.param.fss_ss.noise_temperature = 950
        self.param.fss_ss.antenna_gain = 51
        self.param.fss_ss.antenna_pattern = "OMNI"
        self.param.fss_ss.imt_altitude = 1000
        self.param.fss_ss.imt_lat_deg = -23.5629739
        self.param.fss_ss.imt_long_diff_deg = (-46.6555132-75)
        self.param.fss_ss.channel_model = "FSPL"
        self.param.fss_ss.line_of_sight_prob = 0.01
        self.param.fss_ss.surf_water_vapour_density = 7.5
        self.param.fss_ss.specific_gaseous_att = 0.1
        self.param.fss_ss.time_ratio = 0.5
        self.param.fss_ss.antenna_l_s = -20
        self.param.fss_ss.acs = 0
        self.param.fss_ss.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.fss_ss.EARTH_RADIUS = 6371000

        self.param.fss_es.x = -5000
        self.param.fss_es.y = 0
        self.param.fss_es.location = "FIXED"
        self.param.fss_es.height = 10
        self.param.fss_es.elevation_min = 20
        self.param.fss_es.elevation_max = 20
        self.param.fss_es.azimuth = "0"
        self.param.fss_es.frequency = 10000
        self.param.fss_es.bandwidth = 100
        self.param.fss_es.noise_temperature = 100
        self.param.fss_es.tx_power_density = -60
        self.param.fss_es.antenna_gain = 50
        self.param.fss_es.antenna_pattern = "OMNI"
        self.param.fss_es.channel_model = "FSPL"
        self.param.fss_es.line_of_sight_prob = 1
        self.param.fss_es.acs = 0
        self.param.fss_es.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.fss_es.EARTH_RADIUS = 6371000

        self.param.ras.x = -5000
        self.param.ras.y = 0
        self.param.ras.height = 10
        self.param.ras.elevation = 20
        self.param.ras.azimuth = 0
        self.param.ras.frequency = 10000
        self.param.ras.bandwidth = 100
        self.param.ras.antenna_noise_temperature = 50
        self.param.ras.receiver_noise_temperature = 50
        self.param.ras.antenna_gain = 50
        self.param.ras.antenna_efficiency = 0.7
        self.param.ras.diameter = 10
        self.param.ras.acs = 0
        self.param.ras.antenna_pattern = "OMNI"
        self.param.ras.channel_model = "FSPL"
        self.param.ras.line_of_sight_prob = 1
        self.param.ras.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.ras.EARTH_RADIUS = 6371000
        self.param.ras.SPEED_OF_LIGHT = 299792458


    def test_simulation_2bs_4ue_fss_ss(self):
        self.param.general.system = "FSS_SS"

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        self.assertTrue(self.simulation.co_channel)

        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        self.simulation.ue.x = np.array([20, 70, 110, 170])
        self.simulation.ue.y = np.array([ 0,  0,   0,   0])
        self.simulation.ue.antenna = np.array([AntennaOmni(10), AntennaOmni(11), AntennaOmni(22), AntennaOmni(23)])
        self.simulation.ue.active = np.ones(4, dtype=bool)

        # test connection method
        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        self.simulation.link = {0:[0,1],1:[2,3]}
        self.assertEqual(self.simulation.link, {0: [0,1], 1: [2,3]})

        # We do not test the selection method here because in this specific
        # scenario we do not want to change the order of the UE's

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.fss_ss.channel_model,
                                                                                   self.param, random_number_gen)

        # test coupling loss method
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        path_loss_imt = np.array([[78.47,  89.35,  93.27,  97.05],
                                  [97.55,  94.72,  91.53,  81.99]])
        bs_antenna_gains = np.array([[ 1,  1,  1,  1], [ 2,  2,  2,  2]])
        ue_antenna_gains = np.array([[ 10,  11,  22,  23], [ 10,  11,  22,  23]])
        coupling_loss_imt = path_loss_imt - bs_antenna_gains - ue_antenna_gains \
                            + self.param.imt.bs_ohmic_loss \
                            + self.param.imt.ue_ohmic_loss \
                            + self.param.imt.ue_body_loss

        npt.assert_allclose(self.simulation.coupling_loss_imt,
                            coupling_loss_imt,
                            atol=1e-2)

        # test scheduler and bandwidth allocation
        self.simulation.scheduler()
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)
        npt.assert_allclose(self.simulation.ue.bandwidth, bandwidth_per_ue*np.ones(4), atol=1e-2)

        # there is no power control, so BS's will transmit at maximum power
        self.simulation.power_control()

        tx_power = 10 - 10*math.log10(2)
        npt.assert_allclose(self.simulation.bs.tx_power[0], np.array([tx_power, tx_power]), atol=1e-2)
        npt.assert_allclose(self.simulation.bs.tx_power[1], np.array([tx_power, tx_power]), atol=1e-2)

        # test method that calculates SINR
        self.simulation.calculate_sinr()

        # check UE received power

        rx_power = tx_power - np.concatenate((coupling_loss_imt[0][:2], coupling_loss_imt[1][2:])) 
        npt.assert_allclose(self.simulation.ue.rx_power, rx_power, atol=1e-2)

        # check UE received interference
        rx_interference = tx_power - np.concatenate((coupling_loss_imt[1][:2], coupling_loss_imt[0][2:]))     
        npt.assert_allclose(self.simulation.ue.rx_interference, rx_interference, atol=1e-2)

        # check UE thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9
        npt.assert_allclose(self.simulation.ue.thermal_noise, thermal_noise, atol=1e-2)

        # check UE thermal noise + interference
        total_interference = 10*np.log10(np.power(10, 0.1*rx_interference) + np.power(10, 0.1*thermal_noise))
        npt.assert_allclose(self.simulation.ue.total_interference, total_interference, atol=1e-2)

        # check SNR
        npt.assert_allclose(self.simulation.ue.snr, rx_power - thermal_noise, atol=1e-2)

        # check SINR
        npt.assert_allclose(self.simulation.ue.sinr, rx_power - total_interference, atol=1e-2)

        self.simulation.system = StationFactory.generate_fss_space_station(self.param.fss_ss)
        self.simulation.system.x = np.array([0.01]) # avoids zero-division
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.fss_ss.altitude])

        # test the method that calculates interference from IMT UE to FSS space station
        self.simulation.calculate_external_interference()

        # check coupling loss
        # 4 values because we have 2 BS * 2 beams for each base station. 
        path_loss_imt_system = 203.52
        polarization_loss = 3
        sat_antenna_gain = 51
        bs_antenna_gain = np.array([1, 2])
        coupling_loss_imt_system = path_loss_imt_system - sat_antenna_gain \
                                    - np.array([bs_antenna_gain[0], bs_antenna_gain[0], bs_antenna_gain[1], bs_antenna_gain[1]]) \
                                    + polarization_loss \
                                    + self.param.imt.bs_ohmic_loss
                                    
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        # check interference generated by BS to FSS space station
        interference = tx_power - coupling_loss_imt_system
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
#        self.assertAlmostEqual(self.simulation.system.inr,
#                               np.array([ -147.448 - (-88.821) ]),
#                               delta=.01)
        self.assertAlmostEqual(self.simulation.system.inr,
                               np.array([ rx_interference - thermal_noise ]),
                               delta=.01)

    def test_simulation_2bs_4ue_fss_es(self):
        self.param.general.system = "FSS_ES"

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()


        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        self.simulation.ue.x = np.array([20, 70, 110, 170])
        self.simulation.ue.y = np.array([ 0,  0,   0,   0])
        self.simulation.ue.antenna = np.array([AntennaOmni(10), AntennaOmni(11), AntennaOmni(22), AntennaOmni(23)])
        self.simulation.ue.active = np.ones(4, dtype=bool)

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)

        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        self.simulation.link = {0:[0,1],1:[2,3]}
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        self.simulation.scheduler()
        self.simulation.power_control()
        self.simulation.calculate_sinr()

        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)

        tx_power = 10 - 10*math.log10(2)
        npt.assert_allclose(self.simulation.bs.tx_power[0], np.array([tx_power, tx_power]), atol=1e-2)
        npt.assert_allclose(self.simulation.bs.tx_power[1], np.array([tx_power, tx_power]), atol=1e-2)

        # check UE received power
        rx_power = np.array([tx_power-3-(78.47-1-10)-4-3, tx_power-3-(89.35-1-11)-4-3, tx_power-3-(91.53-2-22)-4-3, tx_power-3-(81.99-2-23)-4-3])
        npt.assert_allclose(self.simulation.ue.rx_power, rx_power, atol=1e-2)

        # check UE received interference
        rx_interference = np.array([tx_power-3-(97.55-2-10)-4-3,  tx_power-3-(94.72-2-11)-4-3, tx_power-3-(93.27-1-22)-4-3, tx_power-3-(97.05-1-23)-4-3])
        npt.assert_allclose(self.simulation.ue.rx_interference, rx_interference, atol=1e-2)

        # check UE thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9
        npt.assert_allclose(self.simulation.ue.thermal_noise, thermal_noise, atol=1e-2)

        # check UE thermal noise + interference
        total_interference = 10*np.log10(np.power(10, 0.1*rx_interference) + np.power(10, 0.1*thermal_noise))
        npt.assert_allclose(self.simulation.ue.total_interference, total_interference, atol=1e-2)

        self.simulation.system = StationFactory.generate_fss_earth_station(self.param.fss_es, random_number_gen)
        self.simulation.system.x = np.array([-2000])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.fss_es.height])

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.fss_es.channel_model,
                                                                                   self.param, random_number_gen)
        # what if FSS ES is the interferer?
        self.simulation.calculate_sinr_ext()

        # check coupling loss between FSS_ES and IMT_UE
        coupling_loss_imt_system = np.array([128.55-50-10,  128.76-50-11,  128.93-50-22,  129.17-50-23])
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        # check interference from FSS_ES to IMT_UE
        system_tx_power = -60 + 10*math.log10(bandwidth_per_ue*1e6) + 30
        ext_interference = system_tx_power - coupling_loss_imt_system
        npt.assert_allclose(self.simulation.ue.ext_interference,
                            ext_interference,
                            atol=1e-2)

        ext_interference_total = 10*np.log10(np.power(10, 0.1*total_interference) \
                                           + np.power(10, 0.1*ext_interference))

        npt.assert_allclose(self.simulation.ue.sinr_ext,
                            rx_power - ext_interference_total,
                            atol=1e-2)

        npt.assert_allclose(self.simulation.ue.inr,
                            ext_interference - thermal_noise,
                            atol=1e-2)

        # what if IMT is interferer?
        self.simulation.calculate_external_interference()

        # check coupling loss from IMT_BS to FSS_ES
        coupling_loss_imt_system = np.array([124.47-50-1,  124.47-50-1,  125.29-50-2,  125.29-50-2]) 
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        interference = tx_power - coupling_loss_imt_system
        rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference)))
        self.assertAlmostEqual(self.simulation.system.rx_interference,
                               rx_interference,
                               delta=.01)

        # check FSS Earth station thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*100*1e3*100*1e6)
        self.assertAlmostEqual(self.simulation.system.thermal_noise,
                               thermal_noise,
                               delta=.01)

        # check INR at FSS Earth station
        self.assertAlmostEqual(self.simulation.system.inr,
                               np.array([ rx_interference - thermal_noise ]),
                               delta=.01)


    def test_simulation_2bs_4ue_ras(self):
        self.param.general.system = "RAS"

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()


        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        self.simulation.ue.x = np.array([20, 70, 110, 170])
        self.simulation.ue.y = np.array([ 0,  0,   0,   0])
        self.simulation.ue.antenna = np.array([AntennaOmni(10), AntennaOmni(11), AntennaOmni(22), AntennaOmni(23)])
        self.simulation.ue.active = np.ones(4, dtype=bool)

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.ras.channel_model,
                                                                                   self.param, random_number_gen)

        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        self.simulation.link = {0:[0,1],1:[2,3]}
        self.simulation.select_ue(random_number_gen)
        self.simulation.link = {0:[0,1],1:[2,3]}
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)
        self.simulation.scheduler()
        self.simulation.power_control()
        self.simulation.calculate_sinr()

        # check UE thermal noise
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9
        npt.assert_allclose(self.simulation.ue.thermal_noise,
                            thermal_noise,
                            atol=1e-2)

        # check SINR
        npt.assert_allclose(self.simulation.ue.sinr,
                            np.array([-70.48 - (-85.49), -80.36 - (-83.19), -70.54 - (-73.15), -60.00 - (-75.82)]),
                            atol=1e-2)

        self.simulation.system = StationFactory.generate_ras_station(self.param.ras)
        self.simulation.system.x = np.array([-2000])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.ras.height])
        self.simulation.system.antenna[0].effective_area = 54.9779

        # Test gain calculation
        gains = self.simulation.calculate_gains(self.simulation.system,self.simulation.bs)
        npt.assert_equal(gains,np.array([[50, 50]]))

        self.simulation.calculate_external_interference()

        polarization_loss = 3
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            np.array([118.47-50-1,  118.47-50-1,  119.29-50-2,  119.29-50-2]) + polarization_loss,
                            atol=1e-2)

        # Test RAS interference
        interference = self.param.imt.bs_conducted_power - 10*np.log10(self.param.imt.ue_k) \
                       - np.array([118.47-50-1,  118.47-50-1,  119.29-50-2,  119.29-50-2]) - polarization_loss
        rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference)))
        self.assertAlmostEqual(self.simulation.system.rx_interference,
                               rx_interference,
                               delta=.01)

        # Test RAS PFD
        pfd = 10*np.log10(10**(rx_interference/10)/54.9779)
        self.assertAlmostEqual(self.simulation.system.pfd,
                            pfd,
                            delta=.01)

        # check RAS station thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*100*1e3*100*1e6)
        self.assertAlmostEqual(self.simulation.system.thermal_noise,
                               thermal_noise,
                               delta=.01)
        # check INR at RAS station
        self.assertAlmostEqual(self.simulation.system.inr,
                               np.array([ rx_interference - (-98.599) ]),
                               delta=.01)

    def test_calculate_bw_weights(self):
        self.param.general.system = "FSS_ES"
        self.simulation = SimulationDownlink(self.param, "")

        bw_imt = 200
        bw_sys = 33.33
        ue_k = 3
        ref_weights = np.array([ 0.5, 0, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 100
        bw_sys = 25
        ue_k = 3
        ref_weights = np.array([ 0.75, 0, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 66.67
        ue_k = 3
        ref_weights = np.array([ 1, 0, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 400
        bw_sys = 200
        ue_k = 3
        ref_weights = np.array([ 1, 0.49, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 133.33
        ue_k = 3
        ref_weights = np.array([ 1, 1, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 150
        ue_k = 3
        ref_weights = np.array([ 1, 1, 0.25])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 150
        bw_sys = 150
        ue_k = 3
        ref_weights = np.array([ 1, 1, 1])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 300
        ue_k = 3
        ref_weights = np.array([ 1, 1, 1])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 50
        ue_k = 2
        ref_weights = np.array([ 0.5, 0])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 100
        bw_sys = 60
        ue_k = 2
        ref_weights = np.array([ 1, 0.2])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 300
        bw_sys = 300
        ue_k = 2
        ref_weights = np.array([ 1, 1])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 100
        bw_sys = 50
        ue_k = 1
        ref_weights = np.array([ 0.5 ])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

        bw_imt = 200
        bw_sys = 180
        ue_k = 1
        ref_weights = np.array([ 0.9])
        weights = self.simulation.calculate_bw_weights(bw_imt, bw_sys, ue_k)
        npt.assert_allclose(ref_weights, weights, atol=1e-2)

if __name__ == '__main__':
    unittest.main()
