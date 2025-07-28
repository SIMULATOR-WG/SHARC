import unittest
import numpy as np
from sharc.satellite.ngso.orbit_model import OrbitModel


class TestOrbitModel(unittest.TestCase):

    def setUp(self):
        self.orbit = OrbitModel(
            Nsp=6,
            Np=8,
            phasing=7.5,
            long_asc=0,
            omega=0,
            delta=52,
            hp=1414,
            ha=1414,
            Mo=0,
        )

    def test_initialization(self):
        self.assertEqual(self.orbit.Nsp, 6)
        self.assertEqual(self.orbit.Np, 8)
        self.assertAlmostEqual(self.orbit.phasing, 7.5)
        self.assertAlmostEqual(self.orbit.long_asc, 0)
        self.assertAlmostEqual(self.orbit.omega, 0)
        self.assertAlmostEqual(self.orbit.delta, 52)
        self.assertAlmostEqual(self.orbit.perigee_alt_km, 1414)
        self.assertAlmostEqual(self.orbit.apogee_alt_km, 1414)
        self.assertAlmostEqual(self.orbit.Mo, 0)

    def test_orbit_parameters(self):
        self.assertAlmostEqual(self.orbit.semi_major_axis, 7792.145)
        self.assertAlmostEqual(self.orbit.eccentricity, 0.0)
        self.assertAlmostEqual(self.orbit.orbital_period_sec, 6845.3519, places=4)
        self.assertAlmostEqual(self.orbit.sat_sep_angle_deg, 60.0)
        self.assertAlmostEqual(self.orbit.orbital_plane_inclination, 45.0)

    def test_mean_anomalies(self):
        self.orbit.get_orbit_positions_time_instant(time_instant_secs=0)
        ma_deg = np.degrees(self.orbit.mean_anomaly.reshape((self.orbit.Np, self.orbit.Nsp)))

        # Check phasing between satellites in the same plane
        r = np.ones((ma_deg.shape[0] - 1, ma_deg.shape[1])) * self.orbit.phasing
        np.testing.assert_array_almost_equal(np.diff(ma_deg, axis=0), r, decimal=4)

        # Check phase between planes
        r = np.ones((ma_deg.shape[0], ma_deg.shape[1] - 1)) * self.orbit.sat_sep_angle_deg
        np.testing.assert_array_almost_equal(np.diff(ma_deg, axis=1), r, decimal=4)

        # Check mean anomaly values with randon time
        self.orbit.get_orbit_positions_random(rng=np.random.RandomState())
        ma_deg = np.unwrap(np.degrees(self.orbit.mean_anomaly.reshape((self.orbit.Np, self.orbit.Nsp))), period=360)

        # Check phasing between satellites in the same plane
        r = np.ones((ma_deg.shape[0] - 1, ma_deg.shape[1])) * self.orbit.phasing
        phasing_diff = np.diff(ma_deg, axis=0)
        phasing_diff[phasing_diff < 0] += 360
        np.testing.assert_array_almost_equal(phasing_diff, r, decimal=4)

        # Check phase between planes
        r = np.ones((ma_deg.shape[0], ma_deg.shape[1] - 1)) * self.orbit.sat_sep_angle_deg
        np.testing.assert_array_almost_equal(np.diff(ma_deg, axis=1), r, decimal=4)


if __name__ == '__main__':
    unittest.main()
