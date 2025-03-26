from pathlib import Path
import unittest
from sharc.parameters.parameters import Parameters
import numpy as np


class ParametersTest(unittest.TestCase):
    """Run Parameter class tests.
    """

    def setUp(self):
        self.parameters = Parameters()
        param_file = Path(__file__).parent.resolve() / \
            'parameters_for_testing.yaml'
        self.parameters.set_file_name(param_file)
        self.parameters.read_params()

    def test_parameters_imt(self):
        """Unit test for ParametersIMT
        """
        self.assertEqual(self.parameters.imt.topology.type, "INDOOR")
        self.assertEqual(
            self.parameters.imt.minimum_separation_distance_bs_ue, 1.3)
        self.assertEqual(self.parameters.imt.interfered_with, False)
        self.assertEqual(self.parameters.imt.frequency, 24360)
        self.assertEqual(self.parameters.imt.bandwidth, 200.5)
        self.assertEqual(self.parameters.imt.rb_bandwidth, 0.181)
        self.assertEqual(self.parameters.imt.spectral_mask, "3GPP E-UTRA")
        self.assertEqual(self.parameters.imt.spurious_emissions, -13.1)
        self.assertEqual(self.parameters.imt.guard_band_ratio, 0.14)

        self.assertEqual(self.parameters.imt.bs.load_probability, 0.2)
        self.assertEqual(self.parameters.imt.bs.conducted_power, 11.1)
        self.assertEqual(self.parameters.imt.bs.height, 6.1)
        self.assertEqual(self.parameters.imt.bs.noise_figure, 10.1)
        self.assertEqual(self.parameters.imt.bs.ohmic_loss, 3.1)
        self.assertEqual(self.parameters.imt.uplink.attenuation_factor, 0.4)
        self.assertEqual(self.parameters.imt.uplink.sinr_min, -10.0)
        self.assertEqual(self.parameters.imt.uplink.sinr_max, 22.0)
        self.assertEqual(self.parameters.imt.ue.k, 3)
        self.assertEqual(self.parameters.imt.ue.k_m, 1)
        self.assertEqual(self.parameters.imt.ue.indoor_percent, 5.0)
        self.assertEqual(
            self.parameters.imt.ue.distribution_type,
            "ANGLE_AND_DISTANCE")
        self.assertEqual(
            self.parameters.imt.ue.distribution_distance,
            "UNIFORM")
        self.assertEqual(self.parameters.imt.ue.azimuth_range, (-70, 90))
        self.assertEqual(self.parameters.imt.ue.tx_power_control, True)
        self.assertEqual(self.parameters.imt.ue.p_o_pusch, -95.0)
        self.assertEqual(self.parameters.imt.ue.alpha, 1.0)
        self.assertEqual(self.parameters.imt.ue.p_cmax, 22.0)
        self.assertEqual(self.parameters.imt.ue.power_dynamic_range, 63.0)
        self.assertEqual(self.parameters.imt.ue.height, 1.5)
        self.assertEqual(self.parameters.imt.ue.noise_figure, 10.0)
        self.assertEqual(self.parameters.imt.ue.ohmic_loss, 3.0)
        self.assertEqual(self.parameters.imt.ue.body_loss, 4.0)
        self.assertEqual(self.parameters.imt.downlink.attenuation_factor, 0.6)
        self.assertEqual(self.parameters.imt.downlink.sinr_min, -10.0)
        self.assertEqual(self.parameters.imt.downlink.sinr_max, 30.0)
        self.assertEqual(self.parameters.imt.channel_model, "FSPL")
        self.assertEqual(self.parameters.imt.los_adjustment_factor, 18.0)
        self.assertEqual(self.parameters.imt.shadowing, False)

        """Test ParametersImtAntenna parameters
        """
        self.assertEqual(
            self.parameters.imt.adjacent_antenna_model,
            "BEAMFORMING")
        self.assertEqual(self.parameters.imt.bs.antenna.normalization, False)
        self.assertEqual(self.parameters.imt.ue.antenna.normalization, False)
        self.assertEqual(self.parameters.imt.bs.antenna.normalization_file,
                         "antenna/beamforming_normalization/bs_norm.npz")
        self.assertEqual(self.parameters.imt.ue.antenna.normalization_file,
                         "antenna/beamforming_normalization/ue_norm.npz")
        self.assertEqual(
            self.parameters.imt.bs.antenna.element_pattern,
            "F1336")
        self.assertEqual(
            self.parameters.imt.ue.antenna.element_pattern,
            "F1336")
        self.assertEqual(
            self.parameters.imt.bs.antenna.minimum_array_gain, -200)
        self.assertEqual(
            self.parameters.imt.ue.antenna.minimum_array_gain, -200)
        self.assertEqual(self.parameters.imt.bs.antenna.downtilt, 6)
        self.assertEqual(self.parameters.imt.bs.antenna.element_max_g, 5)
        self.assertEqual(self.parameters.imt.ue.antenna.element_max_g, 5)
        self.assertEqual(self.parameters.imt.bs.antenna.element_phi_3db, 65)
        self.assertEqual(self.parameters.imt.ue.antenna.element_phi_3db, 90)
        self.assertEqual(self.parameters.imt.bs.antenna.element_theta_3db, 65)
        self.assertEqual(self.parameters.imt.ue.antenna.element_theta_3db, 90)
        self.assertEqual(self.parameters.imt.bs.antenna.n_rows, 8)
        self.assertEqual(self.parameters.imt.ue.antenna.n_rows, 4)
        self.assertEqual(self.parameters.imt.bs.antenna.n_columns, 8)
        self.assertEqual(self.parameters.imt.ue.antenna.n_columns, 4)
        self.assertEqual(
            self.parameters.imt.bs.antenna.element_horiz_spacing, 0.5)
        self.assertEqual(
            self.parameters.imt.ue.antenna.element_horiz_spacing, 0.5)
        self.assertEqual(
            self.parameters.imt.bs.antenna.element_vert_spacing, 0.5)
        self.assertEqual(
            self.parameters.imt.ue.antenna.element_vert_spacing, 0.5)
        self.assertEqual(self.parameters.imt.bs.antenna.element_am, 30)
        self.assertEqual(self.parameters.imt.ue.antenna.element_am, 25)
        self.assertEqual(self.parameters.imt.bs.antenna.element_sla_v, 30)
        self.assertEqual(self.parameters.imt.ue.antenna.element_sla_v, 25)
        self.assertEqual(
            self.parameters.imt.bs.antenna.multiplication_factor, 12)
        self.assertEqual(
            self.parameters.imt.ue.antenna.multiplication_factor, 12)

        """Test ParametersHotspot
        """
        self.assertEqual(self.parameters.imt.topology.hotspot.num_hotspots_per_cell, 1)
        self.assertEqual(self.parameters.imt.topology.hotspot.max_dist_hotspot_ue, 99.9)
        self.assertEqual(self.parameters.imt.topology.hotspot.min_dist_bs_hotspot, 1.2)
        self.assertEqual(self.parameters.imt.topology.hotspot.intersite_distance, 321)
        self.assertEqual(self.parameters.imt.topology.hotspot.num_clusters, 7)
        self.assertEqual(self.parameters.imt.topology.hotspot.wrap_around, True)

        """Test ParametersMacrocell
        """
        self.assertEqual(self.parameters.imt.topology.macrocell.intersite_distance, 543)
        self.assertEqual(self.parameters.imt.topology.macrocell.num_clusters, 7)
        self.assertEqual(self.parameters.imt.topology.macrocell.wrap_around, True)

        """Test ParametersSingleBaseStation
        """
        self.assertEqual(self.parameters.imt.topology.single_bs.cell_radius, 543)
        self.assertEqual(self.parameters.imt.topology.single_bs.intersite_distance,
                         self.parameters.imt.topology.single_bs.cell_radius * 3 / 2)
        self.assertEqual(self.parameters.imt.topology.single_bs.num_clusters, 2)

        """Test ParametersIndoor
        """
        self.assertEqual(self.parameters.imt.topology.indoor.basic_path_loss, "FSPL")
        self.assertEqual(self.parameters.imt.topology.indoor.n_rows, 3)
        self.assertEqual(self.parameters.imt.topology.indoor.n_colums, 2)
        self.assertEqual(self.parameters.imt.topology.indoor.num_imt_buildings, 2)
        self.assertEqual(self.parameters.imt.topology.indoor.street_width, 30.1)
        self.assertEqual(self.parameters.imt.topology.indoor.intersite_distance, 40.1)
        self.assertEqual(self.parameters.imt.topology.indoor.num_cells, 3)
        self.assertEqual(self.parameters.imt.topology.indoor.num_floors, 1)
        self.assertEqual(self.parameters.imt.topology.indoor.ue_indoor_percent, .95)
        self.assertEqual(
            self.parameters.imt.topology.indoor.building_class,
            "THERMALLY_EFFICIENT")

        self.assertEqual(self.parameters.imt.topology.ntn.bs_height, self.parameters.imt.bs.height)
        self.assertEqual(self.parameters.imt.topology.ntn.cell_radius, 123)
        self.assertEqual(self.parameters.imt.topology.ntn.intersite_distance,
                         self.parameters.imt.topology.ntn.cell_radius * np.sqrt(3))
        self.assertEqual(self.parameters.imt.topology.ntn.bs_azimuth, 45)
        self.assertEqual(self.parameters.imt.topology.ntn.bs_elevation, 45)
        self.assertEqual(self.parameters.imt.topology.ntn.num_sectors, 19)

    def test_parameters_fss_ss(self):
        """Test ParametersFssSs
        """
        self.assertEqual(self.parameters.fss_ss.frequency, 43000.0)
        self.assertEqual(self.parameters.fss_ss.bandwidth, 200.0)
        self.assertEqual(self.parameters.fss_ss.tx_power_density, -5.0)
        self.assertEqual(self.parameters.fss_ss.altitude, 35780000.0)
        self.assertEqual(self.parameters.fss_ss.lat_deg, 0.0)
        self.assertEqual(self.parameters.fss_ss.elevation, 270.0)
        self.assertEqual(self.parameters.fss_ss.azimuth, 0.0)
        self.assertEqual(self.parameters.fss_ss.noise_temperature, 950.0)
        self.assertEqual(self.parameters.fss_ss.adjacent_ch_selectivity, 0.0)
        self.assertEqual(self.parameters.fss_ss.antenna_gain, 46.6)
        self.assertEqual(self.parameters.fss_ss.antenna_pattern, "FSS_SS")
        self.assertEqual(self.parameters.fss_ss.earth_station_alt_m, 0.0)
        self.assertEqual(self.parameters.fss_ss.earth_station_lat_deg, 0.0)
        self.assertEqual(
            self.parameters.fss_ss.earth_station_long_diff_deg, 0.0)
        self.assertEqual(self.parameters.fss_ss.season, "SUMMER")
        self.assertEqual(self.parameters.fss_ss.channel_model, "P619")
        self.assertEqual(self.parameters.fss_ss.antenna_l_s, -20.1)
        self.assertEqual(self.parameters.fss_ss.antenna_3_dB, 0.65)

        self.assertEqual(self.parameters.fss_ss.antenna_s1528.antenna_pattern, "ITU-R-S.1528-LEO")
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.slr, 21)
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.antenna_l_s, -20.1)
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.n_side_lobes, 5)
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.l_r, .4)
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.l_t, .4)
        self.assertEqual(self.parameters.fss_ss.antenna_s1528.roll_off, 2)

    def test_parameters_fss_es(self):
        """Test ParametersFssEs
        """
        self.assertEqual(self.parameters.fss_es.location, "UNIFORM_DIST")
        self.assertEqual(self.parameters.fss_es.x, 10000)
        self.assertEqual(self.parameters.fss_es.y, 0)
        self.assertEqual(self.parameters.fss_es.min_dist_to_bs, 10)
        self.assertEqual(self.parameters.fss_es.max_dist_to_bs, 600)
        self.assertEqual(self.parameters.fss_es.height, 6)
        self.assertEqual(self.parameters.fss_es.elevation_min, 48)
        self.assertEqual(self.parameters.fss_es.elevation_max, 80)
        self.assertEqual(self.parameters.fss_es.azimuth, "RANDOM")
        self.assertEqual(self.parameters.fss_es.frequency, 43000)
        self.assertEqual(self.parameters.fss_es.bandwidth, 6)
        self.assertEqual(self.parameters.fss_es.adjacent_ch_selectivity, 0)
        self.assertEqual(self.parameters.fss_es.tx_power_density, -68.3)
        self.assertEqual(self.parameters.fss_es.noise_temperature, 950)
        self.assertEqual(self.parameters.fss_es.antenna_gain, 32)
        self.assertEqual(
            self.parameters.fss_es.antenna_pattern,
            "Modified ITU-R S.465")
        self.assertEqual(self.parameters.fss_es.antenna_envelope_gain, 0)
        self.assertEqual(self.parameters.fss_es.diameter, 1.8)
        self.assertEqual(self.parameters.fss_es.channel_model, "P452")
        self.assertEqual(self.parameters.fss_es.atmospheric_pressure, 935)
        self.assertEqual(self.parameters.fss_es.air_temperature, 300)
        self.assertEqual(self.parameters.fss_es.N0, 352.58)
        self.assertEqual(self.parameters.fss_es.delta_N, 43.127)
        self.assertEqual(self.parameters.fss_es.percentage_p, 0.2)
        self.assertEqual(self.parameters.fss_es.Dct, 70)
        self.assertEqual(self.parameters.fss_es.Dcr, 70)
        self.assertEqual(self.parameters.fss_es.Hte, 20)
        self.assertEqual(self.parameters.fss_es.Hre, 3)
        self.assertEqual(self.parameters.fss_es.tx_lat, -23.55028)
        self.assertEqual(self.parameters.fss_es.rx_lat, -23.17889)
        self.assertEqual(self.parameters.fss_es.polarization, "horizontal")
        self.assertEqual(self.parameters.fss_es.clutter_loss, True)
        self.assertEqual(self.parameters.fss_es.es_position, "ROOFTOP")
        self.assertEqual(self.parameters.fss_es.shadow_enabled, True)
        self.assertEqual(self.parameters.fss_es.building_loss_enabled, True)
        self.assertEqual(self.parameters.fss_es.same_building_enabled, False)
        self.assertEqual(self.parameters.fss_es.diffraction_enabled, False)
        self.assertEqual(
            self.parameters.fss_es.bs_building_entry_loss_type,
            "P2109_FIXED")
        self.assertEqual(
            self.parameters.fss_es.bs_building_entry_loss_prob, 0.75)
        self.assertEqual(
            self.parameters.fss_es.bs_building_entry_loss_value, 35.0)

    def test_parameters_fs(self):
        """Test ParametersFs
        """
        self.assertEqual(self.parameters.fs.x, 1000.0)
        self.assertEqual(self.parameters.fs.y, 0.0)
        self.assertEqual(self.parameters.fs.height, 15.0)
        self.assertEqual(self.parameters.fs.elevation, -10.0)
        self.assertEqual(self.parameters.fs.azimuth, 180.0)
        self.assertEqual(self.parameters.fs.frequency, 27250.0)
        self.assertEqual(self.parameters.fs.bandwidth, 112.0)
        self.assertEqual(self.parameters.fs.noise_temperature, 290.0)
        self.assertEqual(self.parameters.fs.adjacent_ch_selectivity, 20.0)
        self.assertEqual(self.parameters.fs.tx_power_density, -68.3)
        self.assertEqual(self.parameters.fs.antenna_gain, 36.9)
        self.assertEqual(self.parameters.fs.antenna_pattern, "OMNI")
        self.assertEqual(self.parameters.fs.diameter, 0.3)
        self.assertEqual(self.parameters.fs.channel_model, "TerrestrialSimple")

    def test_parameters_haps(self):
        """Test ParametersHaps
        """
        self.assertEqual(self.parameters.haps.frequency, 27251.1)
        self.assertEqual(self.parameters.haps.bandwidth, 200.0)
        self.assertEqual(self.parameters.haps.antenna_gain, 28.1)
        self.assertEqual(self.parameters.haps.eirp_density, 4.4)
        self.assertEqual(self.parameters.haps.tx_power_density,
                         self.parameters.haps.eirp_density - self.parameters.haps.antenna_gain - 60)
        self.assertEqual(self.parameters.haps.altitude, 20001.1)
        self.assertEqual(self.parameters.haps.lat_deg, 0.1)
        self.assertEqual(self.parameters.haps.elevation, 270.0)
        self.assertEqual(self.parameters.haps.azimuth, 0)
        self.assertEqual(self.parameters.haps.antenna_pattern, "OMNI")
        self.assertEqual(self.parameters.haps.earth_station_alt_m, 0.0)
        self.assertEqual(self.parameters.haps.earth_station_lat_deg, 0.0)
        self.assertEqual(self.parameters.haps.earth_station_long_diff_deg, 0.0)
        self.assertEqual(self.parameters.haps.season, "SUMMER")
        self.assertEqual(self.parameters.haps.acs, 30.0)
        self.assertEqual(self.parameters.haps.channel_model, "P619")
        self.assertEqual(self.parameters.haps.antenna_l_n, -25)

    def test_paramters_rns(self):
        """Test ParametersRns
        """
        self.assertEqual(self.parameters.rns.x, 660.1)
        self.assertEqual(self.parameters.rns.y, -370.1)
        self.assertEqual(self.parameters.rns.altitude, 150.1)
        self.assertEqual(self.parameters.rns.frequency, 32000.1)
        self.assertEqual(self.parameters.rns.bandwidth, 60.1)
        self.assertEqual(self.parameters.rns.noise_temperature, 1154.1)
        self.assertEqual(self.parameters.rns.tx_power_density, -70.79)
        self.assertEqual(self.parameters.rns.antenna_gain, 30.1)
        self.assertEqual(self.parameters.rns.antenna_pattern, "ITU-R M.1466")
        self.assertEqual(self.parameters.rns.channel_model, "P619")
        self.assertEqual(self.parameters.rns.season, "SUMMER")
        self.assertEqual(self.parameters.rns.earth_station_alt_m, 0.0)
        self.assertEqual(self.parameters.rns.earth_station_lat_deg, 0.0)
        self.assertEqual(self.parameters.rns.acs, 30.1)

    def test_parametes_ras(self):
        """Test ParametersRas
        """
        self.assertEqual(self.parameters.ras.frequency, 2695)
        self.assertEqual(self.parameters.ras.bandwidth, 10)
        self.assertEqual(self.parameters.ras.noise_temperature, 90)
        self.assertEqual(self.parameters.ras.adjacent_ch_selectivity, 20.1)

        self.assertEqual(self.parameters.ras.geometry.height, 15)
        self.assertEqual(self.parameters.ras.geometry.azimuth.type, "FIXED")
        self.assertEqual(self.parameters.ras.geometry.azimuth.fixed, -90)
        self.assertEqual(self.parameters.ras.geometry.elevation.type, "FIXED")
        self.assertEqual(self.parameters.ras.geometry.elevation.fixed, 45)
        self.assertEqual(self.parameters.ras.antenna.pattern, "OMNI")
        self.assertEqual(self.parameters.ras.antenna.gain, 0.5)
        self.assertEqual(self.parameters.ras.channel_model, "P452")
        self.assertEqual(self.parameters.ras.polarization_loss, 0.0)

    def test_parameters_single_earth_station(self):
        """Test ParametersSingleEarthStation
        """
        self.assertEqual(self.parameters.single_earth_station.frequency, 8250)
        self.assertEqual(self.parameters.single_earth_station.bandwidth, 100)
        self.assertEqual(
            self.parameters.single_earth_station.adjacent_ch_selectivity, 20.0)
        self.assertEqual(
            self.parameters.single_earth_station.tx_power_density, -65.0)
        self.assertEqual(
            self.parameters.single_earth_station.noise_temperature, 300)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.height, 6)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.azimuth.type,
            "FIXED")
        self.assertEqual(
            self.parameters.single_earth_station.geometry.azimuth.fixed, 0)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.azimuth.uniform_dist.min, -180)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.azimuth.uniform_dist.max, 180)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.elevation.type,
            "FIXED")
        self.assertEqual(
            self.parameters.single_earth_station.geometry.elevation.fixed, 60)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.elevation.uniform_dist.min, 30)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.elevation.uniform_dist.max, 65)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.type,
            "CELL")
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.fixed.x, 10)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.fixed.y, 100)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.uniform_dist.min_dist_to_center,
            101)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.uniform_dist.max_dist_to_center,
            102)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.cell.min_dist_to_bs,
            100)
        self.assertEqual(
            self.parameters.single_earth_station.geometry.location.network.min_dist_to_bs,
            150)
        self.assertEqual(self.parameters.single_earth_station.antenna.gain, 28)
        self.assertEqual(
            self.parameters.single_earth_station.antenna.itu_r_f_699.diameter, 1.1)
        self.assertEqual(
            self.parameters.single_earth_station.antenna.itu_r_f_699.frequency,
            self.parameters.single_earth_station.frequency)
        self.assertEqual(
            self.parameters.single_earth_station.antenna.itu_r_f_699.antenna_gain,
            self.parameters.single_earth_station.antenna.gain)

        self.assertEqual(
            self.parameters.single_earth_station.param_p619.earth_station_alt_m, 1200)
        self.assertEqual(
            self.parameters.single_earth_station.param_p619.space_station_alt_m,
            540000)
        self.assertEqual(
            self.parameters.single_earth_station.param_p619.earth_station_lat_deg, 13)
        self.assertEqual(
            self.parameters.single_earth_station.param_p619.earth_station_long_diff_deg,
            10)

        self.assertEqual(
            self.parameters.single_earth_station.param_p452.atmospheric_pressure, 1)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.air_temperature, 2)
        self.assertEqual(self.parameters.single_earth_station.param_p452.N0, 3)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.delta_N, 4)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.percentage_p, 5)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.Dct, 6)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.Dcr, 7)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.Hte, 8)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.Hre, 9)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.tx_lat, 10)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.rx_lat, 11)
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.polarization,
            "horizontal")
        self.assertEqual(
            self.parameters.single_earth_station.param_p452.clutter_loss, True)

        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.es_position,
            "BUILDINGSIDE")
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.shadow_enabled,
            False)
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.building_loss_enabled,
            False)
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.same_building_enabled, True)
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.diffraction_enabled,
            False)
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.bs_building_entry_loss_type,
            "FIXED_VALUE")
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.bs_building_entry_loss_prob,
            0.19)
        self.assertEqual(
            self.parameters.single_earth_station.param_hdfss.bs_building_entry_loss_value,
            47)

        self.parameters.single_earth_station.geometry.azimuth.uniform_dist.max = None
        # this should still not throw, since azimuth is using fixed type
        self.parameters.single_earth_station.validate()

        # now it should throw:
        with self.assertRaises(ValueError) as err_context:
            self.parameters.single_earth_station.geometry.azimuth.type = "UNIFORM_DIST"
            self.parameters.single_earth_station.validate()

        self.assertTrue(
            'azimuth.uniform_dist.max' in str(
                err_context.exception))


if __name__ == '__main__':
    unittest.main()
