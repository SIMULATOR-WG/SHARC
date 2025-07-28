import unittest
import numpy as np
import numpy.testing as npt
import sharc.satellite.utils.sat_utils as sat_utils


class TestSatUtils(unittest.TestCase):

    def setUp(self):
        pass

    def test_ecef2lla(self):
        # Object is over the meridional plane at 1414km of altitude
        sx = 7792137.0
        sy = 0.0
        sz = 0.0
        lat, lng, alt = sat_utils.ecef2lla(sx, sy, sz)
        npt.assert_almost_equal(lat, 0.0, 2)
        npt.assert_almost_equal(lng, 0.0, 2)
        npt.assert_almost_equal(alt, 1414000.0, 1)

        # Object is over the meridional plane at sea level
        sx = 6378137.0
        sy = 0.0
        sz = 0.0
        lat, lng, alt = sat_utils.ecef2lla(sx, sy, sz)
        npt.assert_almost_equal(lat, 0.0, 2)
        npt.assert_almost_equal(lng, 0.0, 2)
        npt.assert_almost_equal(alt, 0.0, 1)

        sx = 7792137.0
        sy = 6378137.0
        sz = 3264751.4
        lat, lng, alt = sat_utils.ecef2lla(sx, sy, sz)
        npt.assert_almost_equal(lat, 18.0316, 3)
        npt.assert_almost_equal(lng, 39.30153, 3)
        npt.assert_almost_equal(alt, 4209582, 1)

        # test the array form
        sx = [7792137.0, 6378137.0, 7792137.0]
        sy = [0.0, 0.0, 6378137.0]
        sz = [0.0, 0.0, 3264751.4]
        lat, lng, alt = sat_utils.ecef2lla(sx, sy, sz)
        npt.assert_almost_equal(lat, [0.0, 0.0, 18.0316], 3)
        npt.assert_almost_equal(lng, [0.0, 0.0, 39.30153], 3)
        npt.assert_almost_equal(alt, [1414000.0, 0.0, 4209582], 1)

    def test_lla2ecef(self):
        # Object is over the meridional plane at 1414km of altitude
        lat = 0.0
        lng = 0.0
        alt = 1414000.0
        sx, sy, sz = sat_utils.lla2ecef(lat, lng, alt)
        npt.assert_almost_equal(sx, 7792137.0, 1)
        npt.assert_almost_equal(sy, 0.0, 1)
        npt.assert_almost_equal(sz, 0.0, 1)

        # Object is over the meridional plane at sea level
        lat = 0.0
        lng = 0.0
        alt = 0.0
        sx, sy, sz = sat_utils.lla2ecef(lat, lng, alt)
        npt.assert_almost_equal(sx, 6378137.0, 1)
        npt.assert_almost_equal(sy, 0.0, 1)
        npt.assert_almost_equal(sz, 0.0, 1)

        lat = 43.0344
        lng = 46.839308
        alt = 141400.0
        sx, sy, sz = sat_utils.lla2ecef(lat, lng, alt)
        npt.assert_almost_equal(sx, 3264751.4, 1)
        npt.assert_almost_equal(sy, 3481390.4, 1)
        npt.assert_almost_equal(sz, 4426792.5, 1)

        # test the array form
        lat = [0.0, 0.0, 43.0344]
        lng = [0.0, 0.0, 46.839308]
        alt = [1414000.0, 0.0, 141400.0]
        sx, sy, sz = sat_utils.lla2ecef(lat, lng, alt)
        npt.assert_almost_equal(sx, [7792137.0, 6378137.0, 3264751.4], 1)
        npt.assert_almost_equal(sy, [0.0, 0.0, 3481390.4], 1)
        npt.assert_almost_equal(sz, [0.0, 0.0, 4426792.5], 1)

    def test_calculate_elev_angle(self):
        earth_station_coords = [
            (0.0, 0.0),
            (-15.0, -42.0),
            (10.0, 20.0),
            (-15.0, -25.0),
        ]

        space_station_coords = [
            (0.0, 0.0),
            (-10.0, -40.0),
            (12.0, 25.0),
            (0.0, -30.0),
        ]

        space_station_alts_km = [
            1414,
            525,
            340,
            35786,
        ]

        expected_elevations = [
            90.0,
            37.48,
            26.67,
            71.45,
        ]

        for i in range(len(expected_elevations)):
            e = sat_utils.calc_elevation(
                earth_station_coords[i][0],
                space_station_coords[i][0],
                earth_station_coords[i][1],
                space_station_coords[i][1],
                space_station_alts_km[i],
            )
            npt.assert_almost_equal(e, expected_elevations[i], 1)


if __name__ == '__main__':
    unittest.main()
