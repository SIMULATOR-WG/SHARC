# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:53:15 2017

@author: Calil

Updated on Sun Jul 29 21:34:07 2024

Request for Issue 33 correction suggestion, adding the unitary tests for LEO type satellite in 600Km and 1200Km heights.

@author: Thiago Ferreira
"""

# Thiago Ferreira - 231025717

# from footprint import Footprint # Old importing for testing and debugging space
# Importing Footprint class from the given module in repository
from sharc.support.footprint import Footprint
import unittest
import numpy as np
import numpy.testing as npt

# Added the matplotlib for plotting the footprint generated in the test
import matplotlib.pyplot as plt


class FootprintTest(unittest.TestCase):
    """
    Unit testing class for Footprint calculations, focusing on different satellite configurations
    such as beam width, elevation angle, and satellite height.
    """

    def setUp(self):
        """
        Set up the test environment by initializing Footprint instances with different parameters.
        """
        # Geostationary (GEO type) satellite height (35786000 m) and different
        # beam widths and elevations
        self.fa1 = Footprint(0.1, bore_lat_deg=0, bore_subsat_long_deg=0.0)
        self.fa2 = Footprint(0.325, bore_lat_deg=0)
        self.fa3 = Footprint(0.325, elevation_deg=20)
        self.fa4 = Footprint(0.325, elevation_deg=30, sat_height=1200000)
        self.fa5 = Footprint(0.325, elevation_deg=30, sat_height=600000)
        self.fa6 = Footprint(
            0.325,
            sat_height=1200000,
            bore_lat_deg=0,
            bore_subsat_long_deg=17.744178387,
        )

        """
        Requested tests for Low Earth Orbit (LEO) Satellite at 1200Km and 600Km added bellow (Issue 33).
        """

        # New tests obtained for the LEO satellite type heights (1200km and
        # 600km)
        # Different satellite heights (1200 km and 600 km)
        self.sat_heights = [1200000, 600000]

    def test_construction(self):
        """
        Test the correct construction of Footprint instances with expected values.
        """
        # Verify properties of fa1
        self.assertEqual(self.fa1.bore_lat_deg, 0)
        self.assertEqual(self.fa1.bore_subsat_long_deg, 0)
        self.assertEqual(self.fa1.beam_width_deg, 0.1)
        self.assertEqual(self.fa1.bore_lat_rad, 0)
        self.assertEqual(self.fa1.bore_subsat_long_rad, 0)
        self.assertEqual(self.fa1.beam_width_rad, np.pi / 1800)
        self.assertEqual(self.fa1.beta, 0)
        self.assertEqual(self.fa1.bore_tilt, 0)

        # Verify properties of fa2
        self.assertEqual(self.fa2.bore_lat_deg, 0)
        self.assertEqual(self.fa2.bore_subsat_long_deg, 0)
        self.assertEqual(self.fa2.bore_lat_rad, 0)
        self.assertEqual(self.fa2.bore_subsat_long_rad, 0)

        # Verify properties of fa3
        self.assertEqual(self.fa3.bore_lat_deg, 0)
        self.assertAlmostEqual(
            self.fa3.bore_subsat_long_deg, 61.84, delta=0.01,
        )

        # Verify properties of fa4
        self.assertEqual(self.fa4.bore_lat_deg, 0)
        self.assertAlmostEqual(
            self.fa4.bore_subsat_long_deg, 13.22, delta=0.01,
        )

        # Verify properties of fa5
        self.assertEqual(self.fa5.bore_lat_deg, 0)
        self.assertAlmostEqual(self.fa5.bore_subsat_long_deg, 7.68, delta=0.01)

        self.assertEqual(self.fa6.sat_height, 1200000)
        self.assertAlmostEqual(self.fa6.elevation_deg, 20, delta=0.01)

    def test_set_elevation(self):
        """
        Test the set_elevation method to ensure it correctly updates the elevation angle.
        """
        self.fa2.set_elevation(20)
        self.assertEqual(self.fa2.bore_lat_deg, 0)
        self.assertAlmostEqual(
            self.fa2.bore_subsat_long_deg, 61.84, delta=0.01,
        )

    def test_calc_footprint(self):
        """
        Test the calc_footprint method to verify the coordinates of the generated footprint polygon.
        """
        fp_long, fp_lat = self.fa1.calc_footprint(4)
        npt.assert_allclose(
            fp_long, np.array(
            [0.0, 0.487, -0.487, 0.0],
            ), atol=1e-2,
        )
        npt.assert_allclose(
            fp_lat, np.array(
            [-0.562, 0.281, 0.281, -0.562],
            ), atol=1e-2,
        )

    def test_calc_area(self):
        """
        Test the calc_area method to verify the calculation of the footprint area.
        """
        a1 = self.fa2.calc_area(1000)
        self.assertAlmostEqual(a1, 130000, delta=130000 * 0.0025)
        a2 = self.fa3.calc_area(1000)
        self.assertAlmostEqual(a2, 486300, delta=486300 * 0.0025)
        a3 = self.fa4.calc_area(1000)
        self.assertAlmostEqual(a3, 810, delta=810 * 0.0025)
        a4 = self.fa5.calc_area(1000)
        self.assertAlmostEqual(a4, 234, delta=234 * 0.0025)

        for height in self.sat_heights:
            beam_deg = 0.325
            footprint = Footprint(
                beam_deg, elevation_deg=90, sat_height=height,
            )
            cone_radius_in_km = height * np.tan(np.deg2rad(beam_deg)) / 1000
            cone_base_area_in_km2 = np.pi * (cone_radius_in_km**2)
            footprint_area_in_km2 = footprint.calc_area(1000)
            self.assertAlmostEqual(
                footprint_area_in_km2,
                cone_base_area_in_km2,
                delta=cone_base_area_in_km2 * 0.01,
            )


if __name__ == '__main__':
    unittest.main()
