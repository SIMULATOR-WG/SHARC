# -*- coding: utf-8 -*-

import unittest
import numpy.testing as npt
import numpy as np
import math

from sharc.antenna.antenna_s1528 import AntennaS1528Taylor
from sharc.antenna.antenna_multiple_transceiver import AntennaMultipleTransceiver
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528
from sharc.station_manager import StationManager


class AntennaAntennaMultipleTransceiverTest(unittest.TestCase):
    def setUp(self):
        param = ParametersAntennaS1528()
        param.antenna_gain = 30
        param.frequency = 2170.0
        param.bandwidth = 5.0
        param.antenna_3_dB_bw = 4.4127

        self.base_antenna = AntennaS1528Taylor(param)

        self.single_antenna = AntennaMultipleTransceiver(
            azimuths=np.array([0.0]),
            elevations=np.array([0.0]),  # so we point at horizon for test
            num_beams=1,
            transceiver_radiation_pattern=self.base_antenna,
        )
        self.double_antenna_pointing_same_way = AntennaMultipleTransceiver(
            azimuths=np.array([0.0, 0.0]),
            elevations=np.array([0.0, 0.0]),  # so we point at horizon for test
            num_beams=2,
            transceiver_radiation_pattern=self.base_antenna,
        )

        altitude = 1000
        cell_radius = 100
        center_dist = 2 * cell_radius
        first_layer_angles = np.linspace(-150, 150, 6)
        sectors7_x = np.concatenate(([0.0], center_dist * np.cos(first_layer_angles)))
        sectors7_y = np.concatenate(([0.0], center_dist * np.sin(first_layer_angles)))

        self.sectors7_azimuth = np.concatenate(([0.0], first_layer_angles))
        self.sectors7_elevation = np.concatenate(([-90.0], np.repeat(np.rad2deg(np.arctan2(center_dist, altitude)) - 90, 6)))

        self.sectors7_antenna = AntennaMultipleTransceiver(
            azimuths=self.sectors7_azimuth,
            elevations=self.sectors7_elevation,
            num_beams=7,
            transceiver_radiation_pattern=self.base_antenna,
        )

    def test_calculate_gain(self):
        """
        We simply compare gains when using antennas separately
        and when using the multiple transceiver
        since P_lin = sum_{i=0}^{n_beams} P_in * loss * g_other_sys * g_transceiver(i) (linear)
        P_lin = P_in * loss * g_other_sys * sum_{i=0}^{n_beams} g_transceiver(i) (linear)
        so we compare the output of the multiple transceiver antenna with
        gain_mult_transceiver = sum_{i=0}^{n_beams} g_transceiver(i)
        """
        phi = np.array([
            0.0, 0.0, 0.0,
        ])
        theta = np.array([
            90.0, 120.0, 150.0,
        ])

        off_axis_angle = np.array([0.0, 30.0, 60.0])
        expected = self.base_antenna.calculate_gain(
            off_axis_angle_vec=off_axis_angle,
            theta_vec=theta,
        )
        actual = self.single_antenna.calculate_gain(
            phi_vec=phi,
            theta_vec=theta,
        )

        self.assertEqual(actual.shape, expected.shape)
        npt.assert_allclose(actual, expected)

        actual_2 = self.double_antenna_pointing_same_way.calculate_gain(
            phi_vec=phi,
            theta_vec=theta,
        )

        self.assertEqual(actual_2.shape, expected.shape)
        # x * 2 (linear) == x + 3.01... (dB)
        npt.assert_allclose(actual_2, expected + 3.01029995664)

        # when calculating gain at (0,0,0), considering antenna at (0,0,altitude),
        # off axis angle will always be elevation + 90
        off_axis_angle = self.sectors7_elevation + 90.
        # and theta will always be 90 - elevation
        theta = 90 - self.sectors7_elevation
        # and phi will always be -azimuth (since antenna is pointing at +az, to reach 0,0 we need -az)
        phi = -self.sectors7_azimuth

        # antenna pointing downwards
        actual7_sec = self.sectors7_antenna.calculate_gain(
            phi_vec=np.array([0.0]),
            theta_vec=np.array([180.0]),
        )

        expected_gains = self.base_antenna.calculate_gain(
                theta_vec=theta,
                off_axis_angle_vec=off_axis_angle,
        )

        expected = np.array([10 * np.log10(np.sum(10**(expected_gains / 10)))])

        self.assertEqual(actual7_sec.shape, expected.shape)
        npt.assert_allclose(actual7_sec, expected)


if __name__ == '__main__':
    unittest.main()
