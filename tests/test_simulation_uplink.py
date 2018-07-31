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
from sharc.parameters.parameters import Parameters
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.station_factory import StationFactory
from sharc.propagation.propagation_factory import PropagationFactory

class SimulationUplinkTest(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()

        self.param.general.imt_link = "UPLINK"
        self.param.general.enable_cochannel = True
        self.param.general.enable_adjacent_channel = False
        self.param.general.overwrite_output = True

        self.param.imt.topology = "SINGLE_BS"
        self.param.imt.num_macrocell_sites = 19
        self.param.imt.num_clusters = 2
        self.param.imt.intersite_distance = 150
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 10000
        self.param.imt.bandwidth = 100
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "ITU 265-E"
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
        self.param.imt.ue_aclr = 20
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
        
        self.param.antenna_imt.normalization = False
        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_normalization_file = None
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
        self.param.antenna_imt.bs_rx_element_phi_deg_3db = 65
        self.param.antenna_imt.bs_rx_element_theta_deg_3db = 65
        self.param.antenna_imt.bs_rx_element_am = 30
        self.param.antenna_imt.bs_rx_element_sla_v = 30
        self.param.antenna_imt.bs_rx_n_rows = 2
        self.param.antenna_imt.bs_rx_n_columns = 2
        self.param.antenna_imt.bs_rx_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_rx_element_vert_spacing = 0.5
        self.param.antenna_imt.bs_downtilt_deg = 10
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_normalization_file = None
        self.param.antenna_imt.ue_tx_element_max_g = 5
        self.param.antenna_imt.ue_tx_element_phi_deg_3db = 65
        self.param.antenna_imt.ue_tx_element_theta_deg_3db = 65
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

        self.param.fss_ss.frequency = 10000
        self.param.fss_ss.bandwidth = 100
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


    def test_simulation_2bs_4ue_ss(self):
        self.param.general.system = "FSS_SS"

        self.simulation = SimulationUplink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0

        self.assertTrue(self.simulation.co_channel)

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
        coupling_loss_imt = np.array([[78.47-1-10,  89.35-1-11,  93.27-1-22,  97.05-1-23],
                                      [97.55-2-10,  94.72-2-11,  91.53-2-22,  81.99-2-23]])
        npt.assert_allclose(self.simulation.coupling_loss_imt,
                            coupling_loss_imt,
                            atol=1e-2)

        # test scheduler and bandwidth allocation
        self.simulation.scheduler()
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)
        npt.assert_allclose(self.simulation.ue.bandwidth, bandwidth_per_ue*np.ones(4), atol=1e-2)

        # there is no power control, so UE's will transmit at maximum power
        self.simulation.power_control()
        tx_power = 20
        npt.assert_allclose(self.simulation.ue.tx_power, tx_power*np.ones(4))

        # test method that calculates SINR
        self.simulation.calculate_sinr()

        # check BS received power
        rx_power = { 0: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[0,0:2]),
                     1: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[1,2:4])}
        npt.assert_allclose(self.simulation.bs.rx_power[0],
                            rx_power[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_power[1],
                            rx_power[1],
                            atol=1e-2)

        # check BS received interference
        rx_interference = { 0: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[0,2:4]),
                            1: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[1,0:2])}

        npt.assert_allclose(self.simulation.bs.rx_interference[0],
                            rx_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_interference[1],
                            rx_interference[1],
                            atol=1e-2)

        # check BS thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 7
        npt.assert_allclose(self.simulation.bs.thermal_noise,
                            thermal_noise,
                            atol=1e-2)

        # check BS thermal noise + interference
        total_interference = { 0: 10*np.log10(np.power(10, 0.1*rx_interference[0]) + np.power(10, 0.1*thermal_noise)),
                               1: 10*np.log10(np.power(10, 0.1*rx_interference[1]) + np.power(10, 0.1*thermal_noise))}
        npt.assert_allclose(self.simulation.bs.total_interference[0],
                            total_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.total_interference[1],
                            total_interference[1],
                            atol=1e-2)

        # check SNR
        npt.assert_allclose(self.simulation.bs.snr[0],
                            rx_power[0] - thermal_noise,
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.snr[1],
                            rx_power[1] - thermal_noise,
                            atol=1e-2)

        # check SINR
        npt.assert_allclose(self.simulation.bs.sinr[0],
                            rx_power[0] - total_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.sinr[1],
                            rx_power[1] - total_interference[1],
                            atol=1e-2)

        self.simulation.system = StationFactory.generate_fss_space_station(self.param.fss_ss)
        self.simulation.system.x = np.array([0])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.fss_ss.altitude])

        # test the method that calculates interference from IMT UE to FSS space station
        self.simulation.calculate_external_interference()

        # check coupling loss
        coupling_loss_imt_system = np.array([203.52-51-10, 203.52-51-11, 203.52-51-22, 203.52-51-23])
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)
        # check interference generated by UE to FSS space station
        interference_ue = tx_power - 3 - 4 - coupling_loss_imt_system
        rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference_ue)))
        self.assertAlmostEqual(self.simulation.system.rx_interference,
                               rx_interference,
                               delta=.01)

        # check FSS space station thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*950*100*1e3*1e6)
        self.assertAlmostEqual(self.simulation.system.thermal_noise,
                               thermal_noise,
                               delta=.01)

        # check INR at FSS space station
        self.assertAlmostEqual(self.simulation.system.inr,
                               rx_interference - thermal_noise,
                               delta=.01)


    def test_simulation_2bs_4ue_es(self):
        self.param.general.system = "FSS_ES"

        self.simulation = SimulationUplink(self.param, "")
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

        self.simulation.connect_ue_to_bs()

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

        self.simulation.scheduler()
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)
        self.simulation.power_control()

        self.simulation.calculate_sinr()

        tx_power = 20

        # check coupling loss IMT
        coupling_loss_imt = np.array([[78.47-1-10,  89.35-1-11,  93.27-1-22,  97.05-1-23],
                                      [97.55-2-10,  94.72-2-11,  91.53-2-22,  81.99-2-23]])
        npt.assert_allclose(self.simulation.coupling_loss_imt,
                            coupling_loss_imt,
                            atol=1e-2)

        # check BS received power
        rx_power = { 0: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[0,0:2]),
                     1: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[1,2:4])}
        npt.assert_allclose(self.simulation.bs.rx_power[0],
                            rx_power[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_power[1],
                            rx_power[1],
                            atol=1e-2)

        # check BS received interference
        rx_interference = { 0: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[0,2:4]),
                            1: np.array([tx_power-3-4-3, tx_power-3-4-3] - coupling_loss_imt[1,0:2])}

        npt.assert_allclose(self.simulation.bs.rx_interference[0],
                            rx_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.rx_interference[1],
                            rx_interference[1],
                            atol=1e-2)

        # check BS thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 7
        npt.assert_allclose(self.simulation.bs.thermal_noise,
                            thermal_noise,
                            atol=1e-2)

        # check BS thermal noise + interference
        total_interference = { 0: 10*np.log10(np.power(10, 0.1*rx_interference[0]) + np.power(10, 0.1*thermal_noise)),
                               1: 10*np.log10(np.power(10, 0.1*rx_interference[1]) + np.power(10, 0.1*thermal_noise))}
        npt.assert_allclose(self.simulation.bs.total_interference[0],
                            total_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.total_interference[1],
                            total_interference[1],
                            atol=1e-2)

        # check SNR
        npt.assert_allclose(self.simulation.bs.snr[0],
                            rx_power[0] - thermal_noise,
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.snr[1],
                            rx_power[1] - thermal_noise,
                            atol=1e-2)

        # check SINR
        npt.assert_allclose(self.simulation.bs.sinr[0],
                            rx_power[0] - total_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.sinr[1],
                            rx_power[1] - total_interference[1],
                            atol=1e-2)

        self.simulation.system = StationFactory.generate_fss_earth_station(self.param.fss_es, random_number_gen)
        self.simulation.system.x = np.array([-2000])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.fss_es.height])

        # what if FSS ES is interferer???
        self.simulation.calculate_sinr_ext()

        # coupling loss FSS_ES <-> IMT BS
        coupling_loss_imt_system = np.array([118.47-50-1,  118.47-50-1,  119.29-50-2,  119.29-50-2])
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        # external interference
        system_tx_power = -60 + 10*math.log10(bandwidth_per_ue*1e6) + 30
        ext_interference = { 0: system_tx_power - coupling_loss_imt_system[0:2] - 3,
                             1: system_tx_power - coupling_loss_imt_system[2:4] - 3}
        npt.assert_allclose(self.simulation.bs.ext_interference[0],
                            ext_interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.ext_interference[1],
                            ext_interference[1],
                            atol=1e-2)

        # SINR with external interference
        interference = { 0: 10*np.log10(np.power(10, 0.1*total_interference[0]) \
                                      + np.power(10, 0.1*ext_interference[0])),
                         1: 10*np.log10(np.power(10, 0.1*total_interference[1]) \
                                      + np.power(10, 0.1*ext_interference[1]))}

        npt.assert_allclose(self.simulation.bs.sinr_ext[0],
                            rx_power[0] - interference[0],
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.sinr_ext[1],
                            rx_power[1] - interference[1],
                            atol=1e-2)

        # INR
        npt.assert_allclose(self.simulation.bs.inr[0],
                            interference[0] - thermal_noise,
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.inr[1],
                            interference[1] - thermal_noise,
                            atol=1e-2)

        # what if IMT is interferer?
        self.simulation.calculate_external_interference()

        # coupling loss
        coupling_loss_imt_system = np.array([118.55-50-10,  118.76-50-11,  118.93-50-22,  119.17-50-23])
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        # interference
        interference = tx_power - 3 - 4 - coupling_loss_imt_system
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

        self.simulation = SimulationUplink(self.param, "")
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

        self.simulation.connect_ue_to_bs()

        # We do not test the selection method here because in this specific
        # scenario we do not want to change the order of the UE's
        #self.simulation.select_ue()

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.fss_ss.channel_model,
                                                                                   self.param, random_number_gen)

        # test coupling loss method
        self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
                                                                                    self.simulation.ue,
                                                                                    self.simulation.propagation_imt)

        self.simulation.scheduler()
        bandwidth_per_ue = math.trunc((1 - 0.1)*100/2)
        self.simulation.power_control()

        self.simulation.calculate_sinr()
        # check BS thermal noise
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 7
        npt.assert_allclose(self.simulation.bs.thermal_noise,
                            thermal_noise,
                            atol=1e-2)

        # check SINR
        npt.assert_allclose(self.simulation.bs.sinr[0],
                            np.array([-57.47 - (-60.27),  -67.35 - (-63.05)]),
                            atol=1e-2)
        npt.assert_allclose(self.simulation.bs.sinr[1],
                            np.array([-57.53 - (-75.41),  -46.99 - (-71.67)]),
                            atol=1e-2)

        # Create system
        self.simulation.system = StationFactory.generate_ras_station(self.param.ras)
        self.simulation.system.x = np.array([-2000])
        self.simulation.system.y = np.array([0])
        self.simulation.system.height = np.array([self.param.ras.height])
        self.simulation.system.antenna[0].effective_area = 54.9779

        # Test gain calculation
        gains = self.simulation.calculate_gains(self.simulation.system,self.simulation.ue)
        npt.assert_equal(gains,np.array([[50, 50, 50, 50]]))

        # Test external interference
        self.simulation.calculate_external_interference()
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            np.array([118.55-50-10,  118.76-50-11,  118.93-50-22,  119.17-50-23]),
                            atol=1e-2)

        # Test RAS PFD
        interference = 20 - np.array([118.55-50-10,  118.76-50-11,  118.93-50-22,  119.17-50-23])- 7
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

    def test_beamforming_gains(self):
        self.param.general.system = "FSS_SS"

        self.simulation = SimulationUplink(self.param, "")
        self.simulation.initialize()

        eps = 1e-2
        random_number_gen = np.random.RandomState()

        # Set scenario
        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology,
                                                                       random_number_gen)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)
        self.simulation.ue.x = np.array([50.000, 43.301, 150.000, 175.000])
        self.simulation.ue.y = np.array([ 0.000, 25.000,   0.000, 43.301])

        # Physical pointing angles
        self.assertEqual(self.simulation.bs.antenna[0].azimuth,0)
        self.assertEqual(self.simulation.bs.antenna[0].elevation,-10)
        self.assertEqual(self.simulation.bs.antenna[1].azimuth,180)
        self.assertEqual(self.simulation.bs.antenna[0].elevation,-10)

        # Change UE pointing
        self.simulation.ue.azimuth = np.array([180, -90, 90, -90])
        self.simulation.ue.elevation = np.array([-30, -15, 15, 30])
        par = self.param.antenna_imt.get_antenna_parameters("UE","TX")
        for i in range(self.simulation.ue.num_stations):
            self.simulation.ue.antenna[i] = AntennaBeamformingImt(par, self.simulation.ue.azimuth[i],
                                                                  self.simulation.ue.elevation[i])
        self.assertEqual(self.simulation.ue.antenna[0].azimuth,180)
        self.assertEqual(self.simulation.ue.antenna[0].elevation,-30)
        self.assertEqual(self.simulation.ue.antenna[1].azimuth,-90)
        self.assertEqual(self.simulation.ue.antenna[1].elevation,-15)
        self.assertEqual(self.simulation.ue.antenna[2].azimuth,90)
        self.assertEqual(self.simulation.ue.antenna[2].elevation,15)
        self.assertEqual(self.simulation.ue.antenna[3].azimuth,-90)
        self.assertEqual(self.simulation.ue.antenna[3].elevation,30)

        # Simulate connection and selection
        self.simulation.connect_ue_to_bs()
        self.assertEqual(self.simulation.link,{0:[0,1],1:[2,3]})

        # Test BS gains
        # Test pointing vector
        phi, theta = self.simulation.bs.get_pointing_vector_to(self.simulation.ue)
        npt.assert_allclose(phi,np.array([[0.0, 30.0, 0.0,    13.898],
                                          [180.0, 170.935, 180.0, 120.0  ]]),atol=eps)
        npt.assert_allclose(theta,np.array([[95.143, 95.143, 91.718, 91.430],
                                            [91.718, 91.624, 95.143, 95.143]]),atol=eps)

        # Add beams by brute force: since the SimulationUplink.select_ue()
        # method shufles the link dictionary, the order of the beams cannot be
        # predicted. Thus, the beams need to be added outside of the function
        self.simulation.ue.active = np.ones(4, dtype=bool)
        self.simulation.bs.antenna[0].add_beam(phi[0,0],theta[0,0])
        self.simulation.bs.antenna[0].add_beam(phi[0,1],theta[0,1])
        self.simulation.bs.antenna[1].add_beam(phi[1,2],theta[1,2])
        self.simulation.bs.antenna[1].add_beam(phi[1,3],theta[1,3])
        self.simulation.ue.antenna[0].add_beam(phi[0,0]-180,180-theta[0,0])
        self.simulation.ue.antenna[1].add_beam(phi[0,1]-180,180-theta[0,1])
        self.simulation.ue.antenna[2].add_beam(phi[1,2]-180,180-theta[1,2])
        self.simulation.ue.antenna[3].add_beam(phi[1,3]-180,180-theta[1,3])
        self.simulation.bs_to_ue_beam_rbs = np.array([0, 1, 0, 1],dtype=int)

        # Test beams pointing
        npt.assert_allclose(self.simulation.bs.antenna[0].beams_list[0],
                           (0.0,-4.857),atol=eps)
        npt.assert_allclose(self.simulation.bs.antenna[0].beams_list[1],
                            (29.92,-3.53),atol=eps)
        npt.assert_allclose(self.simulation.bs.antenna[1].beams_list[0],
                            (0.0,-4.857),atol=eps)
        npt.assert_allclose(self.simulation.bs.antenna[1].beams_list[1],
                            (-59.60,0.10),atol=eps)
        npt.assert_allclose(self.simulation.ue.antenna[0].beams_list[0],
                            (0.0,-35.143),atol=eps)
        npt.assert_allclose(self.simulation.ue.antenna[1].beams_list[0],
                            (-62.04,-12.44),atol=eps)
        npt.assert_allclose(self.simulation.ue.antenna[2].beams_list[0],
                            (-88.66,-4.96),atol=eps)
        npt.assert_allclose(self.simulation.ue.antenna[3].beams_list[0],
                            (32.16,20.71),atol=eps)

        # BS Gain matrix
        ref_gain = np.array([[ 10.954,   8.441,  10.788,   9.460],
                             [ 10.788,   3.365,  10.954,   0.931]])
        gain = self.simulation.calculate_gains(self.simulation.bs,self.simulation.ue)
        npt.assert_allclose(gain,ref_gain,atol=eps)

        # UE Gain matrix
        ref_gain = np.array([[  4.503, -44.198],
                             [ -3.362, -11.206],
                             [-14.812, -14.389],
                             [ -9.726,   3.853]])
        gain = self.simulation.calculate_gains(self.simulation.ue,self.simulation.bs)
        npt.assert_allclose(gain,ref_gain,atol=eps)

    def test_calculate_imt_ul_tput(self):
        self.param.general.system = "FSS_SS"

        self.simulation = SimulationUplink(self.param, "")
        self.simulation.initialize()

        eps = 1e-2

        # Test 1
        snir = np.array([0.0, 1.0, 15.0, -5.0, 100.00, 200.00])
        ref_tput = np.array([ 0.400, 0.470, 2.011, 0.159, 2.927, 2.927])
        tput = self.simulation.calculate_imt_tput(snir,
                                                  self.param.imt.ul_sinr_min,
                                                  self.param.imt.ul_sinr_max,
                                                  self.param.imt.ul_attenuation_factor)
        npt.assert_allclose(tput,ref_tput,atol=eps)

if __name__ == '__main__':
    unittest.main()

