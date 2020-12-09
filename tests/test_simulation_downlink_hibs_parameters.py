"""
Created on Tue Dec 02 11:50:00 2020

@author: Luciano Camilo
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory


class HIBSTestParameters(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()

        self.param.general.imt_link = "DOWNLINK"
        self.param.general.seed = 101
        self.param.general.system = 'RAS'
        self.param.general.enable_cochannel = False
        self.param.general.enable_adjacent_channel = True
        self.param.general.overwrite_output = True

        self.param.imt.topology = "HIBS"
        self.param.imt.wrap_around = False
        self.param.imt.num_clusters = 1
        self.param.imt.intersite_distance = 173205
        self.param.imt.minimum_separation_distance_bs_ue = 0.5
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 2680
        self.param.imt.bandwidth = 20
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "3GPP E-UTRA"
        self.param.imt.spurious_emissions = -13
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.ho_margin = 3
        self.param.imt.bs_load_probability = 0.5
        self.param.imt.num_resource_blocks = 10
        self.param.imt.bs_conducted_power = 37
        self.param.imt.bs_height = 20000
        self.param.imt.bs_acs = 45
        self.param.imt.bs_noise_figure = 5
        self.param.imt.bs_noise_temperature = 290
        self.param.imt.bs_ohmic_loss = 2
        self.param.imt.ul_attenuation_factor = 0.4
        self.param.imt.ul_sinr_min = -10
        self.param.imt.ul_sinr_max = 22
        self.param.imt.ue_k = 3
        self.param.imt.ue_k_m = 1
        self.param.imt.ue_indoor_percent = 0
        self.param.imt.ue_distribution_distance = "RAYLEIGH"
        self.param.imt.ue_distribution_azimuth = "NORMAL"
        self.param.imt.ue_distribution_type = "UNIFORM"
        self.param.imt.ue_tx_power_control = "ON"
        self.param.imt.ue_p_o_pusch = -92
        self.param.imt.ue_alpha = 1
        self.param.imt.ue_p_cmax = 20
        self.param.imt.ue_conducted_power = 23
        self.param.imt.ue_height = 1.5
        self.param.imt.ue_aclr = 20
        self.param.imt.ue_acs = 25
        self.param.imt.ue_noise_figure = 9
        self.param.imt.ue_ohmic_loss = 3
        self.param.imt.ue_body_loss = 4
        self.param.imt.dl_attenuation_factor = 0.8
        self.param.imt.dl_sinr_min = -10
        self.param.imt.dl_sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        self.param.imt.shadowing = True
        self.param.imt.noise_temperature = 290
        self.param.imt.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.imt.EARTH_RADIUS = 6371000

        self.param.antenna_imt.adjacent_antenna_model = "BEAMFORMING"
        self.param.antenna_imt.bs_antenna_type = "BEAMFORMING"
        self.param.antenna_imt.bs_normalization = False
        self.param.antenna_imt.bs_normalization_file = None
        self.param.antenna_imt.bs_element_pattern = "M2101"
        self.param.antenna_imt.bs_minimum_array_gain = -200
        self.param.antenna_imt.bs_element_max_g = 8
        self.param.antenna_imt.bs_element_phi_3db = 65
        self.param.antenna_imt.bs_element_theta_3db = 65
        self.param.antenna_imt.bs_element_am = 30
        self.param.antenna_imt.bs_element_sla_v = 25
        self.param.antenna_imt.bs_n_rows = 7
        self.param.antenna_imt.bs_n_columns = 7
        self.param.antenna_imt.bs_element_horiz_spacing = 0.5
        self.param.antenna_imt.bs_element_vert_spacing = 0.5
        self.param.antenna_imt.bs_multiplication_factor = 12
        self.param.antenna_imt.bs_downtilt = -6

        self.param.antenna_imt.ue_normalization_file = None
        self.param.antenna_imt.ue_normalization = False
        self.param.antenna_imt.ue_element_pattern = "M2101"
        self.param.antenna_imt.ue_minimum_array_gain = -200
        self.param.antenna_imt.ue_element_max_g = -3
        self.param.antenna_imt.ue_element_phi_3db = 90
        self.param.antenna_imt.ue_element_theta_3db = 90
        self.param.antenna_imt.ue_element_am = 30
        self.param.antenna_imt.ue_element_sla_v = 25
        self.param.antenna_imt.ue_n_rows = 2
        self.param.antenna_imt.ue_n_columns = 1
        self.param.antenna_imt.ue_element_horiz_spacing = 0.5
        self.param.antenna_imt.ue_element_vert_spacing = 0.5
        self.param.antenna_imt.ue_multiplication_factor = 12

        self.param.hibs.num_sectors = 7
        self.param.hibs.num_clusters = 1
        self.param.hibs.bs_height = 20000
        self.param.hibs.cell_radius = 100000
        self.param.hibs.intersite_distance = 155884
        self.param.hibs.azimuth3 = '60,180,300'
        self.param.hibs.azimuth7 = '0,0,60,120,180,240,300'
        self.param.hibs.azimuth19 = '0,15,30,45,75,90,105,135,150,165,195,210,225,255,270,285,315,330,345'
        self.param.hibs.elevation3 = '-90,-90,-90'
        self.param.hibs.elevation7 = '-90,-23,-23,-23,-23,-23,-23'
        self.param.hibs.elevation19 = '-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30'
        self.param.hibs.bs_conducted_power = 37
        # For 19 Sectors
        # self.param.hibs.bs_conducted_power = 20.1
        self.param.hibs.bs_backoff_power = 3

        self.param.ras.x = 0
        self.param.ras.y = 0
        self.param.ras.height = 10
        self.param.ras.elevation = 10
        self.param.ras.azimuth = 0
        self.param.ras.frequency = 2695
        self.param.ras.bandwidth = 10
        self.param.ras.antenna_noise_temperature = 12
        self.param.ras.receiver_noise_temperature = 10
        self.param.ras.antenna_gain = 0
        self.param.ras.antenna_efficiency = 0.7
        self.param.ras.diameter = 10
        self.param.ras.acs = 20
        self.param.ras.antenna_pattern = "OMNI"
        self.param.ras.channel_model = "FSPL"
        self.param.ras.line_of_sight_prob = 1
        self.param.ras.BOLTZMANN_CONSTANT = 1.38064852e-23
        self.param.ras.EARTH_RADIUS = 6371000
        self.param.ras.SPEED_OF_LIGHT = 299792458

    def test_power_per_sector_hibs_system1(self):

        """
             Test TX Power of HIBS Base Station (1/3/7 Sector) - System 1

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        if self.param.hibs.num_sectors == 1:
            npt.assert_equal(self.simulation.bs.tx_power, 37)

        elif self.param.hibs.num_sectors == 3:
            for i in range(self.param.hibs.num_sectors):
                npt.assert_equal(self.simulation.bs.tx_power[i], 37)

        elif self.param.hibs.num_sectors == 7:
            npt.assert_equal(self.simulation.bs.tx_power[0], 37)
            for i in range(1, self.param.hibs.num_sectors):
                npt.assert_equal(self.simulation.bs.tx_power[i], 34)

    def test_power_per_sector_hibs_system2(self):

        """
            Test TX Power of HIBS Base Station (1/3/7 Sector) - System 1

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        if self.param.hibs.num_sectors == 19:
            for i in range(self.param.hibs.num_sectors):
                # testing for 20.1dBm
                npt.assert_equal(self.simulation.bs.tx_power[i], 20.1)
                # testing for 11dBm
                # npt.assert_equal(self.simulation.bs.tx_power[i], 11)

    def test_antenna_array_per_sector_hibs_system1(self):

        """
            Test the antenna array configuration of HIBS Base Station (1/3/7 Sector) - System 1

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        num_bs = self.simulation.topology.num_base_stations

        if self.param.hibs.num_sectors == 1:
            npt.assert_equal(self.simulation.bs.antenna[0].n_rows, 2)
            npt.assert_equal(self.simulation.bs.antenna[0].n_cols, 2)

        elif self.param.hibs.num_sectors == 3:
            for i in range(num_bs):
                npt.assert_equal(self.simulation.bs.antenna[i].n_rows, 2)
                npt.assert_equal(self.simulation.bs.antenna[i].n_cols, 2)

        elif self.param.hibs.num_sectors == 7:
            for i in range(num_bs):
                if i == 0:
                    npt.assert_equal(self.simulation.bs.antenna[i].n_rows, 2)
                    npt.assert_equal(self.simulation.bs.antenna[i].n_cols, 2)

                else:
                    npt.assert_equal(self.simulation.bs.antenna[i].n_rows, 4)
                    npt.assert_equal(self.simulation.bs.antenna[i].n_cols, 2)

    def test_antenna_array_per_sector_hibs_system2(self):

        """
            Test the antenna array configuration of HIBS Base Station (19 Sector) - System 2

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        num_bs = self.simulation.topology.num_base_stations

        if self.param.hibs.num_sectors == 19:
            for i in range(num_bs):
                # testing for 7x7 array
                npt.assert_equal(self.simulation.bs.antenna[i].n_rows, 7)
                npt.assert_equal(self.simulation.bs.antenna[i].n_cols, 7)
                # testing for 20 x 20 array
                # npt.assert_equal(self.simulation.bs.antenna[i].n_rows, 7)
                # npt.assert_equal(self.simulation.bs.antenna[i].n_cols, 7)

    def test_azimuth_elevation_hibs_system1(self):

        """
            Test the azimuth/elevation of HIBS Base Station (1/73/7 Sector) - System 1

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        num_bs = self.simulation.topology.num_base_stations

        if self.param.hibs.num_sectors == 1:
            npt.assert_equal(self.simulation.bs.antenna[0].azimuth, 0)
            npt.assert_equal(self.simulation.bs.antenna[0].elevation, -90)

        elif self.param.hibs.num_sectors == 3:
            npt.assert_equal(self.simulation.bs.antenna[0].azimuth, 60)
            npt.assert_equal(self.simulation.bs.antenna[1].azimuth, 180)
            npt.assert_equal(self.simulation.bs.antenna[2].azimuth, 300)

            for i in range(1, num_bs):
                npt.assert_equal(self.simulation.bs.antenna[i].elevation, -90)

        elif self.param.hibs.num_sectors == 7:
            npt.assert_equal(self.simulation.bs.antenna[0].azimuth, 0)
            npt.assert_equal(self.simulation.bs.antenna[1].azimuth, 0)
            npt.assert_equal(self.simulation.bs.antenna[2].azimuth, 60)
            npt.assert_equal(self.simulation.bs.antenna[3].azimuth, 120)
            npt.assert_equal(self.simulation.bs.antenna[4].azimuth, 180)
            npt.assert_equal(self.simulation.bs.antenna[5].azimuth, 240)
            npt.assert_equal(self.simulation.bs.antenna[6].azimuth, 300)

            npt.assert_equal(self.simulation.bs.antenna[0].elevation, -90)
            for i in range(1, num_bs):
                npt.assert_equal(self.simulation.bs.antenna[i].elevation, -23)

    def test_azimuth_elevation_hibs_system2(self):

        """
            Test the azimuth/elevation of HIBS Base Station (19 sector) - System 2

        """

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        num_bs = self.simulation.topology.num_base_stations

        if self.param.hibs.num_sectors == 19:
            npt.assert_equal(self.simulation.bs.antenna[0].azimuth, 0)
            npt.assert_equal(self.simulation.bs.antenna[1].azimuth, 15)
            npt.assert_equal(self.simulation.bs.antenna[2].azimuth, 30)
            npt.assert_equal(self.simulation.bs.antenna[3].azimuth, 45)
            npt.assert_equal(self.simulation.bs.antenna[4].azimuth, 75)
            npt.assert_equal(self.simulation.bs.antenna[5].azimuth, 90)
            npt.assert_equal(self.simulation.bs.antenna[6].azimuth, 105)
            npt.assert_equal(self.simulation.bs.antenna[7].azimuth, 135)
            npt.assert_equal(self.simulation.bs.antenna[8].azimuth, 150)
            npt.assert_equal(self.simulation.bs.antenna[9].azimuth, 165)
            npt.assert_equal(self.simulation.bs.antenna[10].azimuth, 195)
            npt.assert_equal(self.simulation.bs.antenna[11].azimuth, 210)
            npt.assert_equal(self.simulation.bs.antenna[12].azimuth, 225)
            npt.assert_equal(self.simulation.bs.antenna[13].azimuth, 255)
            npt.assert_equal(self.simulation.bs.antenna[14].azimuth, 270)
            npt.assert_equal(self.simulation.bs.antenna[15].azimuth, 285)
            npt.assert_equal(self.simulation.bs.antenna[16].azimuth, 315)
            npt.assert_equal(self.simulation.bs.antenna[17].azimuth, 330)
            npt.assert_equal(self.simulation.bs.antenna[18].azimuth, 345)

            for i in range(1, num_bs):
                npt.assert_equal(self.simulation.bs.antenna[i].elevation, -30)
