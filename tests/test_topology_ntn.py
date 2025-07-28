"""
Created on Tue Dec 28 21:41:00 2024

@author: Vitor Borges
"""

import unittest
import numpy as np
import math
from sharc.topology.topology_ntn import TopologyNTN


class TopologyNTNTest(unittest.TestCase):

    def setUp(self):
        self.bs_height = 1000e3  # meters
        self.bs_azimuth = 45  # degrees
        self.bs_elevation = 45  # degrees
        self.beamwidth = 10
        denominator = math.tan(np.radians(self.beamwidth))
        nominator = math.cos(np.radians(self.bs_elevation))
        self.cell_radius = self.bs_height * denominator / nominator  # meters
        self.intersite_distance = self.cell_radius * np.sqrt(3)  # meters

        self.cos = lambda x: np.cos(np.radians(x))
        self.sin = lambda x: np.sin(np.radians(x))
        self.tan = lambda x: np.tan(np.radians(x))

    def test_single_sector(self):
        topology = TopologyNTN(
            self.intersite_distance,
            self.cell_radius,
            self.bs_height,
            self.bs_azimuth,
            self.bs_elevation,
            num_sectors=1,
        )
        topology.calculate_coordinates()
        expected_x = expected_y = expected_z = [0]
        self.assertListEqual(list(topology.x), expected_x)
        self.assertListEqual(list(topology.y), expected_y)
        self.assertListEqual(list(topology.z), expected_z)
        expected_azimuth = [-135.0]
        for actual_azi, expected_azi in zip(
                topology.azimuth, expected_azimuth,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_azi, expected_azi, places=3)
        expected_elevation = [-45.0]
        for expected_elev, actual_elev in zip(
                topology.elevation, expected_elevation,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_elev, expected_elev, places=3)

    def test_seven_sectors(self):
        topology = TopologyNTN(
            self.intersite_distance,
            self.cell_radius,
            self.bs_height,
            self.bs_azimuth,
            self.bs_elevation,
            num_sectors=7,
        )
        topology.calculate_coordinates()

        d = self.intersite_distance

        # defining expected x, y, z
        expected_x = [0]
        expected_x.extend([d * self.cos(30 + 60 * k) for k in range(6)])
        expected_y = [0]
        expected_y.extend([d * self.sin(30 + 60 * k) for k in range(6)])
        expected_z = [0]
        expected_z.extend([0 for _ in range(6)])

        # testing expected x, y, z
        for actual_x, expec_x in zip(topology.x, expected_x):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_x, expec_x, places=3)
        for actual_y, expec_y in zip(topology.y, expected_y):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_y, expec_y, places=3)
        self.assertListEqual(list(topology.z), expected_z)

        # defining expected azimuth and elevation
        space_x = np.repeat(topology.space_station_x, topology.num_sectors)
        space_y = np.repeat(topology.space_station_y, topology.num_sectors)
        space_z = np.repeat(topology.space_station_z, topology.num_sectors)
        # using expected_x and expected_y so this test is independent of other
        # tests
        expected_azimuth = np.arctan2(
            expected_y - space_y,
            expected_x - space_x,
        ) * 180 / np.pi
        expected_dist_xy = np.sqrt(
            (expected_x - space_x)**2 + (expected_y - space_y)**2,
        )
        # using expected_z so this test is independent of other tests
        expected_elevation = np.arctan2(
            expected_z - space_z,
            expected_dist_xy,
        ) * 180 / np.pi

        # testing expected azimuth and elevation
        for actual_azi, expected_azi in zip(
                topology.azimuth, expected_azimuth,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_azi, expected_azi, places=3)
        for expected_elev, actual_elev in zip(
                topology.elevation, expected_elevation,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_elev, expected_elev, places=3)

    def test_nineteen_sectors(self):
        topology = TopologyNTN(
            self.intersite_distance,
            self.cell_radius,
            self.bs_height,
            self.bs_azimuth,
            self.bs_elevation,
            num_sectors=19,
        )
        topology.calculate_coordinates()

        d = self.intersite_distance

        # defining expected x, y, z
        expected_x = [0]
        expected_y = [0]
        expected_x.extend([d * self.cos(30 + 60 * k) for k in range(6)])
        expected_y.extend([d * self.sin(30 + 60 * k) for k in range(6)])
        for k in range(6):
            # already rotated 30 degrees
            angle = 30 + k * 60
            expected_x.append(2 * d * self.cos(angle))
            expected_y.append(2 * d * self.sin(angle))
            expected_x.append(d * self.cos(angle) + d * self.cos(angle + 60))
            expected_y.append(d * self.sin(angle) + d * self.sin(angle + 60))
        expected_z = [0 for _ in range(19)]

        # testing expected x, y, z
        for actual_x, expec_x in zip(topology.x, expected_x):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_x, expec_x, places=3)
        for actual_y, expec_y in zip(topology.y, expected_y):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_y, expec_y, places=3)
        self.assertListEqual(list(topology.z), expected_z)

        # defining expected azimuth and elevation
        space_x = np.repeat(topology.space_station_x, topology.num_sectors)
        space_y = np.repeat(topology.space_station_y, topology.num_sectors)
        space_z = np.repeat(topology.space_station_z, topology.num_sectors)
        # using expected_x and expected_y so this test is independent of other
        # tests
        expected_azimuth = np.arctan2(
            expected_y - space_y,
            expected_x - space_x,
        ) * 180 / np.pi
        expected_distance_xy = np.sqrt(
            (expected_x - space_x)**2 + (expected_y - space_y)**2,
        )
        # using expected_z so this test is independent of other tests
        expected_elevation = np.arctan2(
            expected_z - space_z,
            expected_distance_xy,
        ) * 180 / np.pi

        # testing expected azimuth and elevation
        for actual_azi, expected_azi in zip(
                topology.azimuth, expected_azimuth,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_azi, expected_azi, places=3)
        for expected_elev, actual_elev in zip(
                topology.elevation, expected_elevation,
        ):
            # cannot check asserListEqual because of float precision
            self.assertAlmostEqual(actual_elev, expected_elev, places=3)


if __name__ == '__main__':
    unittest.main()
