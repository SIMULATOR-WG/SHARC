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
from sharc.parameters.imt.parameters_imt_topology import ParametersImtTopology
from sharc.parameters.imt.parameters_single_bs import ParametersSingleBS


class SimulationDownlinkHapsTest(unittest.TestCase):

    def setUp(self):
        self.param = Parameters()

        self.param.general.imt_link = "DOWNLINK"
        self.param.general.seed = 101
        self.param.general.enable_cochannel = True
        self.param.general.enable_adjacent_channel = False
        self.param.general.overwrite_output = True

        self.param.imt.topology = ParametersImtTopology(
            type="SINGLE_BS",
            single_bs=ParametersSingleBS(
                num_clusters=2,
                intersite_distance=150,
                cell_radius=2 * 150 / 3,
            ),
        )
        self.param.imt.minimum_separation_distance_bs_ue = 10
        self.param.imt.interfered_with = False
        self.param.imt.frequency = 10000.0
        self.param.imt.bandwidth = 100
        self.param.imt.rb_bandwidth = 0.180
        self.param.imt.spectral_mask = "IMT-2020"
        self.param.imt.spurious_emissions = -13
        self.param.imt.guard_band_ratio = 0.1
        self.param.imt.ho_margin = 3
        self.param.imt.bs.load_probability = 1

        self.param.imt.bs.conducted_power = 10
        self.param.imt.bs.height = 6
        self.param.imt.bs.acs = 30
        self.param.imt.bs.noise_figure = 7
        self.param.imt.bs.ohmic_loss = 3
        self.param.imt.uplink.attenuation_factor = 0.4
        self.param.imt.uplink.sinr_min = -10
        self.param.imt.uplink.sinr_max = 22
        self.param.imt.ue.k = 2
        self.param.imt.ue.k_m = 1
        self.param.imt.ue.indoor_percent = 0
        self.param.imt.ue.distribution_distance = "RAYLEIGH"
        self.param.imt.ue.distribution_azimuth = "UNIFORM"
        self.param.imt.ue.distribution_type = "ANGLE_AND_DISTANCE"
        self.param.imt.ue.tx_power_control = "OFF"
        self.param.imt.ue.p_o_pusch = -95
        self.param.imt.ue.alpha = 0.8
        self.param.imt.ue.p_cmax = 20
        self.param.imt.ue.conducted_power = 10
        self.param.imt.ue.height = 1.5
        self.param.imt.ue.aclr = 20
        self.param.imt.ue.acs = 25
        self.param.imt.ue.noise_figure = 9
        self.param.imt.ue.ohmic_loss = 3
        self.param.imt.ue.body_loss = 4
        self.param.imt.downlink.attenuation_factor = 0.6
        self.param.imt.downlink.sinr_min = -10
        self.param.imt.downlink.sinr_max = 30
        self.param.imt.channel_model = "FSPL"
        # probability of line-of-sight (not for FSPL)
        self.param.imt.line_of_sight_prob = 0.75
        self.param.imt.shadowing = False
        self.param.imt.noise_temperature = 290

        self.param.imt.bs.antenna.type = "ARRAY"
        self.param.imt.bs.antenna.array.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.imt.ue.antenna.array.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.imt.bs.antenna.array.normalization = False
        self.param.imt.bs.antenna.array.normalization_file = None
        self.param.imt.bs.antenna.array.element_pattern = "M2101"
        self.param.imt.bs.antenna.array.minimum_array_gain = -200
        self.param.imt.bs.antenna.array.element_max_g = 10
        self.param.imt.bs.antenna.array.element_phi_3db = 80
        self.param.imt.bs.antenna.array.element_theta_3db = 80
        self.param.imt.bs.antenna.array.element_am = 25
        self.param.imt.bs.antenna.array.element_sla_v = 25
        self.param.imt.bs.antenna.array.n_rows = 16
        self.param.imt.bs.antenna.array.n_columns = 16
        self.param.imt.bs.antenna.array.element_horiz_spacing = 1
        self.param.imt.bs.antenna.array.element_vert_spacing = 1
        self.param.imt.bs.antenna.array.multiplication_factor = 12
        self.param.imt.bs.antenna.array.downtilt = 10

        self.param.imt.ue.antenna.type = "ARRAY"
        self.param.imt.ue.antenna.array.normalization_file = None
        self.param.imt.ue.antenna.array.normalization = False
        self.param.imt.ue.antenna.array.element_pattern = "M2101"
        self.param.imt.ue.antenna.array.minimum_array_gain = -200
        self.param.imt.ue.antenna.array.element_max_g = 5
        self.param.imt.ue.antenna.array.element_phi_3db = 65
        self.param.imt.ue.antenna.array.element_theta_3db = 65
        self.param.imt.ue.antenna.array.element_am = 30
        self.param.imt.ue.antenna.array.element_sla_v = 30
        self.param.imt.ue.antenna.array.n_rows = 2
        self.param.imt.ue.antenna.array.n_columns = 1
        self.param.imt.ue.antenna.array.element_horiz_spacing = 0.5
        self.param.imt.ue.antenna.array.element_vert_spacing = 0.5
        self.param.imt.ue.antenna.array.multiplication_factor = 12

        self.param.haps.frequency = 10000
        self.param.haps.bandwidth = 200
        self.param.haps.altitude = 20000
        self.param.haps.lat_deg = 0
        self.param.haps.elevation = 270
        self.param.haps.azimuth = 0
        self.param.haps.eirp_density = 4.4
        self.param.haps.antenna_gain = 28
        self.param.haps.tx_power_density = self.param.haps.eirp_density - \
            self.param.haps.antenna_gain - 60
        self.param.haps.antenna_pattern = "OMNI"
        self.param.haps.imt_altitude = 0
        self.param.haps.imt_lat_deg = 0
        self.param.haps.imt_long_diff_deg = 0
        self.param.haps.season = "SUMMER"
        self.param.haps.channel_model = "FSPL"
        self.param.haps.antenna_l_n = -25

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

        self.simulation.bs = StationFactory.generate_imt_base_stations(
            self.param.imt,
            self.param.imt.bs.antenna.array,
            self.simulation.topology,
            random_number_gen,
        )
        self.simulation.bs.antenna = np.array([AntennaOmni(1), AntennaOmni(2)])
        self.simulation.bs.active = np.ones(2, dtype=bool)

        self.simulation.ue = StationFactory.generate_imt_ue(
            self.param.imt,
            self.param.imt.ue.antenna.array,
            self.simulation.topology,
            random_number_gen,
        )
        self.simulation.ue.x = np.array([20, 70, 110, 170])
        self.simulation.ue.y = np.array([0, 0, 0, 0])
        self.simulation.ue.antenna = np.array(
            [AntennaOmni(10), AntennaOmni(11), AntennaOmni(22), AntennaOmni(23)],
        )
        self.simulation.ue.active = np.ones(4, dtype=bool)

        self.simulation.propagation_imt = PropagationFactory.create_propagation(
            self.param.imt.channel_model,
            self.param,
            self.simulation.param_system,
            random_number_gen,
        )
        self.simulation.propagation_system = PropagationFactory.create_propagation(
            self.param.haps.channel_model,
            self.param,
            self.simulation.param_system,
            random_number_gen,
        )
        self.simulation.connect_ue_to_bs()
        self.simulation.select_ue(random_number_gen)
        self.simulation.link = {0: [0, 1], 1: [2, 3]}
        self.simulation.coupling_loss_imt = \
            self.simulation.calculate_intra_imt_coupling_loss(
                self.simulation.ue,
                self.simulation.bs,
            )
        self.simulation.scheduler()
        self.simulation.power_control()
        self.simulation.calculate_sinr()

        bandwidth_per_ue = math.trunc((1 - 0.1) * 100 / 2)
        thermal_noise = 10 * \
            np.log10(1.38064852e-23 * 290 * bandwidth_per_ue * 1e3 * 1e6) + 9

        tx_power = 10 - 10 * math.log10(2)
        npt.assert_allclose(
            self.simulation.bs.tx_power[0], np.array(
            [tx_power, tx_power],
            ), atol=1e-2,
        )
        npt.assert_allclose(
            self.simulation.bs.tx_power[1], np.array(
            [tx_power, tx_power],
            ), atol=1e-2,
        )

        # check UE received power
        rx_power = np.array([
            tx_power - 3 - (78.68 - 1 - 10) - 4 - 3,
            tx_power - 3 - (89.37 - 1 - 11) - 4 - 3,
            tx_power - 3 - (91.54 - 2 - 22) - 4 - 3,
            tx_power - 3 - (82.09 - 2 - 23) - 4 - 3,
        ])
        npt.assert_allclose(self.simulation.ue.rx_power, rx_power, atol=1e-2)

        # check UE received interference
        rx_interference = np.array([
            tx_power - 3 - (97.55 - 2 - 10) - 4 - 3,
            tx_power - 3 - (94.73 - 2 - 11) - 4 - 3,
            tx_power - 3 - (93.28 - 1 - 22) - 4 - 3,
            tx_power - 3 - (97.07 - 1 - 23) - 4 - 3,
        ])
        npt.assert_allclose(
            self.simulation.ue.rx_interference,
            rx_interference,
            atol=1e-2,
        )

        # check UE thermal noise
        ue_noise_fig = 9
        self.simulation.ue.noise_figure = np.ones(self.simulation.ue.num_stations) * ue_noise_fig
        thermal_noise = 10 * \
            np.log10(1.38064852e-23 * 290 * bandwidth_per_ue * 1e3 * 1e6) + ue_noise_fig
        # the simulator adds noise figure to total noise
        npt.assert_allclose(
            self.simulation.ue.thermal_noise,
            thermal_noise,
            atol=1e-2,
        )

        # check UE thermal noise + interference
        total_interference = 10 * \
            np.log10(np.power(10, 0.1 * rx_interference) + np.power(10, 0.1 * thermal_noise))
        npt.assert_allclose(
            self.simulation.ue.total_interference,
            total_interference,
            atol=1e-2,
        )

        self.simulation.system = StationFactory.generate_haps(
            self.param.haps, 0, random_number_gen,
        )

        # now we evaluate interference from HAPS to IMT UE
        self.simulation.calculate_sinr_ext()

        # check coupling loss between FSS_ES and IMT_UE
        coupling_loss_imt_system = np.array(
            [148.47 - 28 - 10, 148.47 - 28 - 11, 148.47 - 28 - 22, 148.47 - 28 - 23],
        ).reshape((-1, 1))

        npt.assert_allclose(
            self.simulation.coupling_loss_imt_system,
            coupling_loss_imt_system,
            atol=1e-2,
        )

        system_tx_power = (4.4 - 28 - 60) + 10 * \
            np.log10(bandwidth_per_ue * 1e6) + 30

        ext_interference = (system_tx_power - coupling_loss_imt_system).flatten()
        npt.assert_allclose(
            self.simulation.ue.ext_interference,
            ext_interference,
            atol=1e-2,
        )

        ext_interference_total = 10 * np.log10(
            np.power(10, 0.1 * total_interference) +
            np.power(10, 0.1 * ext_interference),
        )

        npt.assert_allclose(
            self.simulation.ue.sinr_ext,
            rx_power - ext_interference_total,
            atol=1e-2,
        )

        npt.assert_allclose(
            self.simulation.ue.inr,
            ext_interference - thermal_noise,
            atol=1e-2,
        )


if __name__ == '__main__':
    unittest.main()
