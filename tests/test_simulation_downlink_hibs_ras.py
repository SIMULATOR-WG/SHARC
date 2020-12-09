"""
Created on Tue Dec 03 11:50:00 2020

@author: Luciano Camilo
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.propagation.propagation_factory import PropagationFactory
from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory


class SimulationTestHIBSRAS(unittest.TestCase):

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

        self.param.antenna_imt.adjacent_antenna_model = "SINGLE_ELEMENT"
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
        self.param.antenna_imt.bs_n_rows = 2
        self.param.antenna_imt.bs_n_columns = 2
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

        self.param.hibs.num_sectors = 1
        self.param.hibs.num_clusters = 1
        self.param.hibs.bs_height = 20000
        self.param.hibs.cell_radius = 100000
        self.param.hibs.intersite_distance = 173205
        self.param.hibs.azimuth3 = '60,180,300'
        self.param.hibs.azimuth7 = '0,0,60,120,180,240,300'
        self.param.hibs.azimuth19 = '0,15,30,45,75,90,105,135,150,165,195,210,225,255,270,285,315,330,345'
        self.param.hibs.elevation3 = '-90,-90,-90'
        self.param.hibs.elevation7 = '-90,-23,-23,-23,-23,-23,-23'
        self.param.hibs.elevation19 = '-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30'
        self.param.hibs.bs_conducted_power = 37
        self.param.hibs.bs_backoff_power = 3

        self.param.ras.x = 0
        self.param.ras.y = 0
        self.param.ras.height = 0
        self.param.ras.elevation = 90
        self.param.ras.azimuth = 0
        self.param.ras.frequency = 2695
        self.param.ras.bandwidth = 10
        self.param.ras.antenna_noise_temperature = 12
        self.param.ras.receiver_noise_temperature = 10
        self.param.ras.antenna_gain = 0
        self.param.ras.antenna_efficiency = 1
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

    def test_simulation_hibs_1sector_ras(self):
        """
            Test the interference generated by one HIBS (1 sector) station to
            one RAS base station
        """
        self.param.general.system = "RAS"

        self.simulation = SimulationDownlink(self.param, "")
        self.simulation.initialize()

        self.simulation.bs_power_gain = 0
        self.simulation.ue_power_gain = 0

        random_number_gen = np.random.RandomState()

        self.simulation.bs = StationFactory.generate_imt_base_stations(self.param.imt,
                                                                       self.param.hibs,
                                                                       self.param.antenna_imt,
                                                                       self.simulation.topology, random_number_gen)

        self.simulation.ue = StationFactory.generate_imt_ue(self.param.imt,
                                                            self.param.antenna_imt,
                                                            self.simulation.topology,
                                                            random_number_gen)

        self.simulation.propagation_imt = PropagationFactory.create_propagation(self.param.imt.channel_model,
                                                                                self.param, random_number_gen)
        self.simulation.propagation_system = PropagationFactory.create_propagation(self.param.ras.channel_model,
                                                                                   self.param, random_number_gen)

        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        #self.simulation.coupling_loss_imt = self.simulation.calculate_coupling_loss(self.simulation.bs,
        #                                                                            self.simulation.ue,
        #                                                                            self.simulation.propagation_imt)
        self.simulation.scheduler()
        self.simulation.power_control()


        self.simulation.system = StationFactory.generate_ras_station(self.param.ras)

        # now we evaluate interference from HIBS to RAS
        self.simulation.calculate_sinr()
        self.simulation.calculate_external_interference()

        print(f"Path Loss {self.simulation.imt_system_path_loss}")
        print(f"Perda por acoplamento {self.simulation.coupling_loss_imt_system_adjacent}")
        print(f"Valor de RX Interference {self.simulation.system.rx_interference}")
        print(f"Valor de Ruido Térmico {self.simulation.system.thermal_noise}")

        print(f"RX Interference - Ruido termico {self.simulation.system.rx_interference-self.simulation.system.thermal_noise}")
        print(f"Valor de System INR {self.simulation.system.inr}")
        print(f"Valor de System PFD {self.simulation.system.pfd}")
