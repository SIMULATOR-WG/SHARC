# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:59:07 2017

@author: edgar
"""

import unittest

from sharc.support.enumerations import StationType
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.station import Station

class StationTest(unittest.TestCase):

    def setUp(self):
        #Array parameters
        self.param = ParametersAntennaImt()

        self.param.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.bs_normalization = False
        self.param.bs_normalization_file = None
        self.param.bs_element_pattern = "M2101"
        self.param.bs_minimum_array_gain = -200
        self.param.bs_element_max_g = 10
        self.param.bs_element_phi_3db = 65
        self.param.bs_element_theta_3db = 75
        self.param.bs_element_am = 35
        self.param.bs_element_sla_v = 25
        self.param.bs_n_rows = 8
        self.param.bs_n_columns = 8
        self.param.bs_element_horiz_spacing = 0.5
        self.param.bs_element_vert_spacing = 0.5
        self.param.bs_multiplication_factor = 12
        self.param.bs_downtilt = 0

        self.param.ue_element_pattern = "M2101"
        self.param.ue_normalization = False
        self.param.ue_normalization_file = None
        self.param.ue_minimum_array_gain = -200
        self.param.ue_element_max_g = 10
        self.param.ue_element_phi_3db = 75
        self.param.ue_element_theta_3db = 65
        self.param.ue_element_am = 25
        self.param.ue_element_sla_v = 35
        self.param.ue_n_rows = 2
        self.param.ue_n_columns = 2
        self.param.ue_element_horiz_spacing = 0.5
        self.param.ue_element_vert_spacing = 0.5
        self.param.ue_multiplication_factor = 12

        self.station = Station()
        self.station.id = 1
        self.station.station_type = StationType.IMT_BS
        self.station.x = 10
        self.station.y = 15
        self.station.height = 6
        self.station.tx_power = 20
        self.station.rx_power = -3
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.station.antenna = AntennaBeamformingImt(par,300,-10)

        self.station2 = Station()
        self.station2.id = 1
        self.station2.station_type = StationType.IMT_UE
        self.station2.x = 10
        self.station2.y = 15
        self.station2.height = 6
        self.station2.tx_power = 17
        self.station2.rx_power = 9
        par = self.param.get_antenna_parameters(StationType.IMT_UE)
        self.station2.antenna = AntennaBeamformingImt(par,270,2)

        self.station3 = Station()
        self.station3.id = 2
        self.station3.station_type = StationType.FSS_SS
        self.station3.x = 10
        self.station3.y = 15
        self.station3.height = 6
        self.station3.tx_power = 20
        self.station3.rx_power = -3
        self.station3.antenna = AntennaOmni(50)

    def test_id(self):
        self.assertEqual(self.station.id, 1)

    def test_station_type(self):
        self.assertEqual(self.station.station_type,StationType.IMT_BS)

    def test_x(self):
        self.assertEqual(self.station.x, 10)

    def test_y(self):
        self.assertEqual(self.station.y, 15)

    def test_height(self):
        self.assertEqual(self.station.height, 6)

    def test_tx_power(self):
        self.assertEqual(self.station.tx_power, 20)

    def test_rx_power(self):
        self.assertEqual(self.station.rx_power, -3)

    def test_antenna(self):
        self.assertEqual(self.station.antenna.azimuth, 300)
        self.assertEqual(self.station.antenna.elevation, -10)
        self.assertEqual(self.station.antenna.n_rows, 8)
        self.assertEqual(self.station.antenna.n_cols, 8)
        self.assertEqual(self.station.antenna.dh, 0.5)
        self.assertEqual(self.station.antenna.dv, 0.5)
        self.assertEqual(self.station.antenna.element.g_max, 10)
        self.assertEqual(self.station.antenna.element.phi_3db, 65)
        self.assertEqual(self.station.antenna.element.theta_3db, 75)
        self.assertEqual(self.station.antenna.element.am, 35)
        self.assertEqual(self.station.antenna.element.sla_v, 25)

    def test_eq(self):
        self.assertTrue(self.station == self.station2)
        # changing id, x, y, or height should change the result
        self.station.x = 11
        self.assertFalse(self.station == self.station2)
        #
        self.assertFalse(self.station == self.station3)
        self.assertFalse(self.station2 == self.station3)

    def test_ne(self):
        self.assertFalse(self.station != self.station2)
        # changing id, x, y, or height should change the result
        self.station.height = 9
        self.assertTrue(self.station != self.station2)
        #
        self.assertTrue(self.station != self.station3)
        self.assertTrue(self.station2 != self.station3)


if __name__ == '__main__':
    unittest.main()
