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

class SimulationDownlinkHapsTest(unittest.TestCase):

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

        self.param.antenna_imt.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.antenna_imt.normalization = False
        self.param.antenna_imt.bs_normalization_file = None
        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_minimum_array_gain = -200
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
        self.param.antenna_imt.ue_normalization_file = None
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_minimum_array_gain = -200
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



    def test_simulation_2bs_4ue_1haps(self):
        """
        Test the interference generated by one HAPS (airbone) station to
        one IMT base station
        """
        self.param.general.system = "HAPS"

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
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.haps.channel_model,
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
        thermal_noise = 10*np.log10(1.38064852e-23*290*bandwidth_per_ue*1e3*1e6) + 9

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

        self.simulation.system = StationFactory.generate_haps(self.param.haps, 0, random_number_gen)

        # now we evaluate interference from HAPS to IMT UE
        self.simulation.calculate_sinr_ext()

        # check coupling loss between FSS_ES and IMT_UE
        coupling_loss_imt_system = np.array([148.47-28-10,  148.47-28-11,  148.47-28-22,  148.47-28-23])
        npt.assert_allclose(self.simulation.coupling_loss_imt_system,
                            coupling_loss_imt_system,
                            atol=1e-2)

        system_tx_power = (4.4 - 28 - 60) + 10*math.log10(bandwidth_per_ue*1e6) + 30

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





if __name__ == '__main__':
    unittest.main()
