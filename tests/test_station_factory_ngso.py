import unittest
from sharc.parameters.parameters_mss_d2d import ParametersOrbit, ParametersMssD2d
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.topology.topology_single_base_station import TopologySingleBaseStation
from sharc.support.enumerations import StationType
from sharc.station_factory import StationFactory
from sharc.station_manager import StationManager
from sharc.support.sharc_geom import GeometryConverter, lla2ecef

import numpy as np
import numpy.testing as npt


class StationFactoryNgsoTest(unittest.TestCase):
    def setUp(self):
        # Adding multiple shells to this constellation
        # Creating orbital parameters for the first orbit
        orbit_1 = ParametersOrbit(
            n_planes=20,                  # Number of orbital planes
            sats_per_plane=32,            # Satellites per plane
            phasing_deg=3.9,              # Phasing angle in degrees
            long_asc_deg=18.0,            # Longitude of ascending node
            inclination_deg=54.5,         # Orbital inclination in degrees
            perigee_alt_km=525.0,         # Perigee altitude in kilometers
            apogee_alt_km=525.0,           # Apogee altitude in kilometers
        )

        # Creating orbital parameters for the second orbit
        orbit_2 = ParametersOrbit(
            n_planes=12,                  # Number of orbital planes
            sats_per_plane=20,            # Satellites per plane
            phasing_deg=2.0,              # Phasing angle in degrees
            long_asc_deg=30.0,            # Longitude of ascending node
            inclination_deg=26.0,         # Orbital inclination in degrees
            perigee_alt_km=580.0,         # Perigee altitude in kilometers
            apogee_alt_km=580.0,           # Apogee altitude in kilometers
        )

        # Creating an NGSO constellation and adding the defined orbits
        self.lat = -15.7801
        self.long = -47.9292
        self.alt = 1200

        self.geoconvert = GeometryConverter()
        self.geoconvert.set_reference(
            -15.7801,
            -47.9292,
            1200,
        )
        self.param = ParametersMssD2d(
            name="Acme-Star-1",                         # Name of the constellation
            antenna_pattern="ITU-R-S.1528-Taylor",     # Antenna type
            orbits=[orbit_1, orbit_2],                   # List of orbital parameters
            num_sectors=1,
        )
        self.param.antenna_s1528.frequency = 43000.0
        self.param.antenna_s1528.bandwidth = 500.0
        self.param.antenna_s1528.antenna_gain = 46.6

        # Creating an IMT topology
        imt_topology = TopologySingleBaseStation(
            cell_radius=500,
            num_clusters=2,
        )

        # random number generator
        self.seed = 42
        rng = np.random.RandomState(seed=self.seed)

        self.ngso_manager = StationFactory.generate_mss_d2d(self.param, rng, self.geoconvert)

    def test_ngso_manager(self):
        self.assertEqual(self.ngso_manager.station_type, StationType.MSS_D2D)
        self.assertEqual(self.ngso_manager.num_stations, 20 * 32 + 12 * 20)
        self.assertEqual(self.ngso_manager.x.shape, (20 * 32 + 12 * 20,))
        self.assertEqual(self.ngso_manager.y.shape, (20 * 32 + 12 * 20,))
        self.assertEqual(self.ngso_manager.height.shape, (20 * 32 + 12 * 20,))

    def test_satellite_antenna_pointing(self):
        # by default, satellites should always point to nadir (earth center)

        # Test: check if azimuth is pointing towards correct direction
        # y > 0 <=> azimuth < 0
        # y < 0 <=> azimuth > 0
        npt.assert_array_equal(np.sign(self.ngso_manager.azimuth), -np.sign(self.ngso_manager.y))

        # Test: check if center of earth is 0deg off axis, and that its distance to satellite is correct
        earth_center = StationManager(1)
        earth_center.x = np.array([0.])
        earth_center.y = np.array([0.])
        x, y, z = lla2ecef(self.lat, self.long, self.alt)
        earth_center.z = np.array([
            -np.sqrt(
                x * x + y * y + z * z,
            ),
        ])

        self.assertNotAlmostEqual(earth_center.z[0], 0.)

        off_axis_angle = self.ngso_manager.get_off_axis_angle(earth_center)
        distance_to_center_of_earth = self.ngso_manager.get_3d_distance_to(earth_center)
        distance_to_center_of_earth_should_eq = np.sqrt(
            self.ngso_manager.x ** 2 +
            self.ngso_manager.y ** 2 +
            (
                    np.sqrt(
                        x * x + y * y + z * z,
                    ) + self.ngso_manager.z
            ) ** 2,
        )

        npt.assert_allclose(off_axis_angle, 0.0, atol=1e-05)

        npt.assert_allclose(
            distance_to_center_of_earth.flatten(),
            distance_to_center_of_earth_should_eq,
            atol=1e-05,
        )

    def test_satellite_coordinate_reversing(self):
        # by default, satellites should always point to nadir (earth center)
        rng = np.random.RandomState(seed=self.seed)

        ngso_original_coord = StationFactory.generate_mss_d2d(self.param, rng, self.geoconvert)
        self.geoconvert.revert_station_2d_to_3d(ngso_original_coord)
        # Test: check if azimuth is pointing towards correct direction
        # y > 0 <=> azimuth < 0
        # y < 0 <=> azimuth > 0
        npt.assert_array_equal(np.sign(ngso_original_coord.azimuth), -np.sign(ngso_original_coord.y))

        # Test: check if center of earth is 0deg off axis
        earth_center = StationManager(1)
        earth_center.x = np.array([0.])
        earth_center.y = np.array([0.])
        earth_center.z = np.array([0.])

        off_axis_angle = ngso_original_coord.get_off_axis_angle(earth_center)

        npt.assert_allclose(off_axis_angle, 0.0, atol=1e-05)

        self.geoconvert.convert_station_3d_to_2d(ngso_original_coord)

        npt.assert_allclose(self.ngso_manager.x, ngso_original_coord.x, atol=1e-500)
        npt.assert_allclose(self.ngso_manager.y, ngso_original_coord.y, atol=1e-500)
        npt.assert_allclose(self.ngso_manager.z, ngso_original_coord.z, atol=1e-500)
        npt.assert_allclose(self.ngso_manager.height, ngso_original_coord.height, atol=1e-500)
        npt.assert_allclose(self.ngso_manager.azimuth, ngso_original_coord.azimuth, atol=1e-500)
        npt.assert_allclose(self.ngso_manager.elevation, ngso_original_coord.elevation, atol=1e-500)


if __name__ == '__main__':
    unittest.main()
