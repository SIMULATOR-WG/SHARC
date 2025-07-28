import unittest
import numpy as np
import numpy.testing as npt
from sharc.topology.topology_imt_mss_dc import TopologyImtMssDc
from sharc.parameters.imt.parameters_imt_mss_dc import ParametersImtMssDc
from sharc.station_manager import StationManager
from sharc.parameters.parameters_orbit import ParametersOrbit
from sharc.support.sharc_geom import GeometryConverter, lla2ecef


class TestTopologyImtMssDc(unittest.TestCase):
    """
    Unit tests for the TopologyImtMssDc class.
    This test suite includes the following tests:
    - test_initialization: Verifies the initialization of the IMT MSS-DC topology.
    - test_calculate_coordinates: Tests the calculation of coordinates for the topology.
    - test_visible_satellites: Checks the visibility of satellites based on elevation angle.
    Classes:
        TestTopologyImtMssDc: Contains unit tests for the TopologyImtMssDc class.
    Methods:
        setUp: Sets up the test environment, including parameters and geometry converter.
        test_initialization: Tests the initialization of the IMT MSS-DC topology.
        test_calculate_coordinates: Tests the calculation of coordinates for the topology.
        test_visible_satellites: Tests the visibility of satellites based on elevation angle.
    """

    def setUp(self):
        # Define the parameters for the IMT MSS-DC topology
        orbit = ParametersOrbit(
            n_planes=20,
            sats_per_plane=32,
            phasing_deg=3.9,
            long_asc_deg=18.0,
            inclination_deg=54.5,
            perigee_alt_km=525,
            apogee_alt_km=525,
        )
        self.params = ParametersImtMssDc(
            beam_radius=36516.0,
            num_beams=19,
            orbits=[orbit],
        )
        self.params.sat_is_active_if.conditions = ["MINIMUM_ELEVATION_FROM_ES"]
        self.params.sat_is_active_if.minimum_elevation_from_es = 5.0

        # Define the geometry converter
        self.geometry_converter = GeometryConverter()
        self.geometry_converter.set_reference(-15.0, -42.0, 1200)

        # Define the Earth center coordinates
        self.earth_center_x = np.array([0.])
        self.earth_center_y = np.array([0.])
        x, y, z = lla2ecef(
            self.geometry_converter.ref_lat,
            self.geometry_converter.ref_long,
            self.geometry_converter.ref_alt,
        )
        self.earth_center_z = np.array([-np.sqrt(x * x + y * y + z * z)])

        # Instantiate the IMT MSS-DC topology
        self.imt_mss_dc_topology = TopologyImtMssDc(self.params, self.geometry_converter)

    def test_initialization(self):
        self.assertTrue(self.imt_mss_dc_topology.is_space_station)
        self.assertEqual(self.imt_mss_dc_topology.num_sectors, self.params.num_beams)
        self.assertEqual(len(self.imt_mss_dc_topology.orbits), len(self.params.orbits))

    def test_calculate_coordinates(self):
        self.imt_mss_dc_topology.calculate_coordinates()
        center_beam_idxs = np.arange(self.imt_mss_dc_topology.num_base_stations // self.imt_mss_dc_topology.num_sectors) *\
            self.imt_mss_dc_topology.num_sectors
        self.assertIsNotNone(self.imt_mss_dc_topology.space_station_x)
        self.assertIsNotNone(self.imt_mss_dc_topology.space_station_y)
        self.assertIsNotNone(self.imt_mss_dc_topology.space_station_z)
        self.assertIsNotNone(self.imt_mss_dc_topology.elevation)
        self.assertIsNotNone(self.imt_mss_dc_topology.azimuth)
        self.assertEqual(len(self.imt_mss_dc_topology.space_station_x), self.imt_mss_dc_topology.num_base_stations)
        self.assertEqual(len(self.imt_mss_dc_topology.space_station_y), self.imt_mss_dc_topology.num_base_stations)
        self.assertEqual(len(self.imt_mss_dc_topology.space_station_z), self.imt_mss_dc_topology.num_base_stations)

        # Test: check if azimuth is pointing towards correct direction
        # y > 0 <=> azimuth < 0
        # y < 0 <=> azimuth > 0
        npt.assert_array_equal(
            np.sign(self.imt_mss_dc_topology.azimuth[center_beam_idxs]),
            -np.sign(self.imt_mss_dc_topology.space_station_y[center_beam_idxs]),
        )

        # Test: check if the altitude is calculated correctly
        rx = self.imt_mss_dc_topology.space_station_x - self.earth_center_x
        ry = self.imt_mss_dc_topology.space_station_y - self.earth_center_y
        rz = self.imt_mss_dc_topology.space_station_z - self.earth_center_z
        r = np.sqrt(rx**2 + ry**2 + rz**2)
        expected_alt_km = (r - np.abs(self.earth_center_z)) / 1e3
        npt.assert_array_almost_equal(expected_alt_km, self.params.orbits[0].apogee_alt_km, decimal=0)

        # by default, satellites should always point to nadir (earth center)
        ref_earth_center = StationManager(1)
        ref_earth_center.x = self.earth_center_x
        ref_earth_center.y = self.earth_center_y
        ref_earth_center.z = self.earth_center_z

        ref_space_stations = StationManager(self.imt_mss_dc_topology.num_base_stations)
        ref_space_stations.x = self.imt_mss_dc_topology.space_station_x
        ref_space_stations.y = self.imt_mss_dc_topology.space_station_y
        ref_space_stations.z = self.imt_mss_dc_topology.space_station_z

        phi, theta = ref_space_stations.get_pointing_vector_to(ref_earth_center)
        npt.assert_array_almost_equal(
            np.squeeze(phi[center_beam_idxs]), self.imt_mss_dc_topology.azimuth[center_beam_idxs],
            decimal=3,
        )
        npt.assert_array_almost_equal(
            np.squeeze(theta[center_beam_idxs]), 90 - self.imt_mss_dc_topology.elevation[center_beam_idxs],
            decimal=3,
        )

    def test_visible_satellites(self):
        self.imt_mss_dc_topology.calculate_coordinates(random_number_gen=np.random.RandomState(8))
        min_elevation_angle = 5.0

        # calculate the elevation angles with respect to the x-y plane
        xy_plane_elevations = np.degrees(
            np.arctan2(
                self.imt_mss_dc_topology.space_station_z,
                np.sqrt(self.imt_mss_dc_topology.space_station_x**2 + self.imt_mss_dc_topology.space_station_y**2),
            ),
        )
        # Add a tolerance to the elevation angle because of the Earth oblateness
        expected_atol = 2e-2
        npt.assert_array_less(min_elevation_angle - expected_atol, xy_plane_elevations)


if __name__ == '__main__':
    unittest.main()
