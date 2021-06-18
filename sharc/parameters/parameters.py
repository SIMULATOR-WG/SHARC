# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 19:35:52 2017

@author: edgar
@modified: Luciano Camilo Tue Jun 18 07:20:00 2021
"""

import configparser

from sharc.parameters.parameters_general import ParametersGeneral
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_hotspot import ParametersHotspot
from sharc.parameters.parameters_indoor import ParametersIndoor
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_eess_passive import ParametersEessPassive
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_ss_mleo import ParametersSsMLeo
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_arns import ParametersArns
from sharc.parameters.parameters_ras import ParametersRas
from sharc.parameters.parameters_hibs import ParametersHibs
from sharc.parameters.parameters_imt_base_station import ParametersImtBaseStation
from sharc.parameters.parameters_antenna_imt_base_station import ParametersAntennaImtBaseStation


class Parameters(object):
    """
    Reads parameters from input file.
    """

    def __init__(self):
        self.file_name = None

        self.general = ParametersGeneral()
        self.imt = ParametersImt()
        self.antenna_imt = ParametersAntennaImt()
        self.hotspot = ParametersHotspot()
        self.indoor = ParametersIndoor()
        self.eess_passive = ParametersEessPassive()
        self.fs = ParametersFs()
        self.ss_mleo = ParametersSsMLeo()
        self.fss_ss = ParametersFssSs()
        self.fss_es = ParametersFssEs()
        self.haps = ParametersHaps()
        self.rns = ParametersRns()
        self.arns = ParametersArns()
        self.ras = ParametersRas()
        self.hibs = ParametersHibs()
        self.imt_bs = ParametersImtBaseStation()
        self.antenna_imt_bs = ParametersAntennaImtBaseStation()

    def set_file_name(self, file_name: str):
        self.file_name = file_name

    def read_params(self):
        config = configparser.ConfigParser()
        config.read(self.file_name)

        #######################################################################
        # GENERAL
        #######################################################################
        self.general.num_snapshots           = config.getint("GENERAL", "num_snapshots")
        self.general.imt_link                = config.get("GENERAL", "imt_link")
        self.general.system                  = config.get("GENERAL", "system")
        self.general.enable_cochannel        = config.getboolean("GENERAL", "enable_cochannel")
        self.general.enable_adjacent_channel = config.getboolean("GENERAL", "enable_adjacent_channel")
        self.general.seed                    = config.getint("GENERAL", "seed")
        self.general.overwrite_output        = config.getboolean("GENERAL", "overwrite_output")

        #######################################################################
        # IMT
        #######################################################################
        self.imt.topology                 = config.get("IMT", "topology")
        self.imt.wrap_around              = config.getboolean("IMT", "wrap_around")
        self.imt.num_clusters             = config.getint("IMT", "num_clusters")
        self.imt.intersite_distance       = config.getfloat("IMT", "intersite_distance")
        self.imt.minimum_separation_distance_bs_ue = config.getfloat("IMT", "minimum_separation_distance_bs_ue")
        self.imt.interfered_with          = config.getboolean("IMT", "interfered_with")
        self.imt.frequency                = config.getfloat("IMT", "frequency")
        self.imt.bandwidth                = config.getfloat("IMT", "bandwidth")
        self.imt.rb_bandwidth             = config.getfloat("IMT", "rb_bandwidth")
        self.imt.spectral_mask            = config.get("IMT", "spectral_mask")
        self.imt.spurious_emissions       = config.getfloat("IMT", "spurious_emissions")
        self.imt.guard_band_ratio         = config.getfloat("IMT", "guard_band_ratio")
        self.imt.bs_load_probability      = config.getfloat("IMT", "bs_load_probability")
        self.imt.bs_conducted_power       = config.getfloat("IMT", "bs_conducted_power")
        self.imt.bs_height                = config.getfloat("IMT", "bs_height")
        self.imt.bs_noise_figure          = config.getfloat("IMT", "bs_noise_figure")
        self.imt.bs_noise_temperature     = config.getfloat("IMT", "bs_noise_temperature")
        self.imt.bs_ohmic_loss            = config.getfloat("IMT", "bs_ohmic_loss")
        self.imt.ul_attenuation_factor    = config.getfloat("IMT", "ul_attenuation_factor")
        self.imt.ul_sinr_min              = config.getfloat("IMT", "ul_sinr_min")
        self.imt.ul_sinr_max              = config.getfloat("IMT", "ul_sinr_max")
        self.imt.ue_k                     = config.getint("IMT", "ue_k")
        self.imt.ue_k_m                   = config.getint("IMT", "ue_k_m")
        self.imt.ue_indoor_percent        = config.getfloat("IMT", "ue_indoor_percent")
        self.imt.ue_distribution_type     = config.get("IMT", "ue_distribution_type")
        self.imt.ue_distribution_distance = config.get("IMT", "ue_distribution_distance")
        self.imt.ue_distribution_azimuth  = config.get("IMT", "ue_distribution_azimuth")
        self.imt.ue_tx_power_control      = config.get("IMT", "ue_tx_power_control")
        self.imt.ue_p_o_pusch             = config.getfloat("IMT", "ue_p_o_pusch")
        self.imt.ue_alpha                 = config.getfloat("IMT", "ue_alpha")
        self.imt.ue_p_cmax                = config.getfloat("IMT", "ue_p_cmax")
        self.imt.ue_power_dynamic_range   = config.getfloat("IMT", "ue_power_dynamic_range")
        self.imt.ue_height                = config.getfloat("IMT", "ue_height")
        self.imt.ue_noise_figure          = config.getfloat("IMT", "ue_noise_figure")
        self.imt.ue_ohmic_loss            = config.getfloat("IMT", "ue_ohmic_loss")
        self.imt.ue_body_loss             = config.getfloat("IMT", "ue_body_loss")
        self.imt.dl_attenuation_factor    = config.getfloat("IMT", "dl_attenuation_factor")
        self.imt.dl_sinr_min              = config.getfloat("IMT", "dl_sinr_min")
        self.imt.dl_sinr_max              = config.getfloat("IMT", "dl_sinr_max")
        self.imt.channel_model            = config.get("IMT", "channel_model")
        self.imt.los_adjustment_factor    = config.getfloat("IMT", "los_adjustment_factor")
        self.imt.shadowing                = config.getboolean("IMT", "shadowing")
        self.imt.noise_temperature        = config.getfloat("IMT", "noise_temperature")
        self.imt.BOLTZMANN_CONSTANT       = config.getfloat("IMT", "BOLTZMANN_CONSTANT")
        self.imt.EARTH_RADIUS             = config.getfloat("IMT", "EARTH_RADIUS")

        #######################################################################
        # IMT ANTENNA
        #######################################################################
        self.antenna_imt.bs_antenna_type            = config.get("IMT_ANTENNA", "bs_antenna_type")
        self.antenna_imt.adjacent_antenna_model     = config.get("IMT_ANTENNA", "adjacent_antenna_model")
        self.antenna_imt.bs_normalization           = config.getboolean("IMT_ANTENNA", "bs_normalization")
        self.antenna_imt.ue_normalization           = config.getboolean("IMT_ANTENNA", "ue_normalization")
        self.antenna_imt.bs_normalization_file      = config.get("IMT_ANTENNA", "bs_normalization_file")
        self.antenna_imt.ue_normalization_file      = config.get("IMT_ANTENNA", "ue_normalization_file")
        self.antenna_imt.bs_element_pattern         = config.get("IMT_ANTENNA", "bs_element_pattern")
        self.antenna_imt.ue_element_pattern         = config.get("IMT_ANTENNA", "ue_element_pattern")

        self.antenna_imt.bs_element_max_g           = config.getfloat("IMT_ANTENNA", "bs_element_max_g")
        self.antenna_imt.bs_element_phi_3db         = config.getfloat("IMT_ANTENNA", "bs_element_phi_3db")
        self.antenna_imt.bs_element_theta_3db       = config.getfloat("IMT_ANTENNA", "bs_element_theta_3db")
        self.antenna_imt.bs_element_am              = config.getfloat("IMT_ANTENNA", "bs_element_am")
        self.antenna_imt.bs_element_sla_v           = config.getfloat("IMT_ANTENNA", "bs_element_sla_v")
        self.antenna_imt.bs_n_rows                  = config.getfloat("IMT_ANTENNA", "bs_n_rows")
        self.antenna_imt.bs_n_columns               = config.getfloat("IMT_ANTENNA", "bs_n_columns")
        self.antenna_imt.bs_element_horiz_spacing   = config.getfloat("IMT_ANTENNA", "bs_element_horiz_spacing")
        self.antenna_imt.bs_element_vert_spacing    = config.getfloat("IMT_ANTENNA", "bs_element_vert_spacing")
        self.antenna_imt.bs_multiplication_factor   = config.getfloat("IMT_ANTENNA", "bs_multiplication_factor")
        self.antenna_imt.bs_minimum_array_gain      = config.getfloat("IMT_ANTENNA", "bs_minimum_array_gain")

        self.antenna_imt.ue_element_max_g           = config.getfloat("IMT_ANTENNA", "ue_element_max_g")
        self.antenna_imt.ue_element_phi_3db         = config.getfloat("IMT_ANTENNA", "ue_element_phi_3db")
        self.antenna_imt.ue_element_theta_3db       = config.getfloat("IMT_ANTENNA", "ue_element_theta_3db")
        self.antenna_imt.ue_element_am              = config.getfloat("IMT_ANTENNA", "ue_element_am")
        self.antenna_imt.ue_element_sla_v           = config.getfloat("IMT_ANTENNA", "ue_element_sla_v")
        self.antenna_imt.ue_n_rows                  = config.getfloat("IMT_ANTENNA", "ue_n_rows")
        self.antenna_imt.ue_n_columns               = config.getfloat("IMT_ANTENNA", "ue_n_columns")
        self.antenna_imt.ue_element_horiz_spacing   = config.getfloat("IMT_ANTENNA", "ue_element_horiz_spacing")
        self.antenna_imt.ue_element_vert_spacing    = config.getfloat("IMT_ANTENNA", "ue_element_vert_spacing")
        self.antenna_imt.ue_multiplication_factor   = config.getfloat("IMT_ANTENNA", "ue_multiplication_factor")
        self.antenna_imt.ue_minimum_array_gain      = config.getfloat("IMT_ANTENNA", "ue_minimum_array_gain")
        self.antenna_imt.bs_downtilt            = config.getfloat("IMT_ANTENNA", "bs_downtilt")

        self.antenna_imt.bf_enable                      = config.get("IMT_ANTENNA", "bf_enable")

        #######################################################################
        # HIBS
        #######################################################################
        self.hibs.num_sectors        = config.getint("HIBS", "num_sectors")
        self.hibs.num_clusters       = config.getint("HIBS", "num_clusters")
        self.hibs.bs_height          = config.getint("HIBS", "bs_height")
        self.hibs.intersite_distance = config.getfloat("HIBS", "intersite_distance")
        self.hibs.cell_radius        = config.getfloat("HIBS", "cell_radius")
        self.hibs.azimuth3           = config.get ('HIBS', "azimuth3")
        self.hibs.azimuth7           = config.get ('HIBS', "azimuth7")
        self.hibs.azimuth19          = config.get ('HIBS', "azimuth19")
        self.hibs.elevation3         = config.get ('HIBS', "elevation3")
        self.hibs.elevation7         = config.get ('HIBS', "elevation7")
        self.hibs.elevation19        = config.get ('HIBS', "elevation19")
        self.hibs.bs_conducted_power = config.getfloat ('HIBS', 'bs_conducted_power')
        self.hibs.bs_backoff_power   = config.getfloat ('HIBS', 'bs_backoff_power')

        #######################################################################
        # HOTSPOT
        #######################################################################
        self.hotspot.num_hotspots_per_cell = config.getint("HOTSPOT", "num_hotspots_per_cell")
        self.hotspot.max_dist_hotspot_ue   = config.getfloat("HOTSPOT", "max_dist_hotspot_ue")
        self.hotspot.min_dist_bs_hotspot   = config.getfloat("HOTSPOT", "min_dist_bs_hotspot")

        #######################################################################
        # INDOOR
        #######################################################################
        self.indoor.basic_path_loss        = config.get("INDOOR", "basic_path_loss")
        self.indoor.n_rows                 = config.getint("INDOOR", "n_rows")
        self.indoor.n_colums               = config.getint("INDOOR", "n_colums")
        self.indoor.num_imt_buildings      = config.get("INDOOR", "num_imt_buildings")
        self.indoor.street_width           = config.getint("INDOOR", "street_width")
        self.indoor.intersite_distance     = config.getfloat("INDOOR", "intersite_distance")
        self.indoor.num_cells              = config.getint("INDOOR", "num_cells")
        self.indoor.num_floors             = config.getint("INDOOR", "num_floors")
        self.indoor.ue_indoor_percent      = config.getfloat("INDOOR", "ue_indoor_percent")
        self.indoor.building_class         = config.get("INDOOR", "building_class")

        #######################################################################
        # SS Medium/Low earth orbit space station
        #######################################################################
        self.ss_mleo.frequency = config.getfloat("SS_MLEO", "frequency")
        self.ss_mleo.bandwidth = config.getfloat("SS_MLEO", "bandwidth")
        self.ss_mleo.tx_power_density = config.getfloat("FSS_SS", "tx_power_density")
        self.ss_mleo.altitude = config.getfloat("SS_MLEO", "altitude")
        self.ss_mleo.lat_deg = config.getfloat("SS_MLEO", "lat_deg")
        self.ss_mleo.elevation = config.getfloat("SS_MLEO", "elevation")
        self.ss_mleo.azimuth = config.getfloat("SS_MLEO", "azimuth")
        self.ss_mleo.noise_temperature = config.getfloat("SS_MLEO", "noise_temperature")
        self.ss_mleo.adjacent_ch_selectivity = config.getfloat("SS_MLEO", "adjacent_ch_selectivity")
        self.ss_mleo.antenna_gain = config.getfloat("SS_MLEO", "antenna_gain")
        self.ss_mleo.antenna_pattern = config.get("SS_MLEO", "antenna_pattern")
        self.ss_mleo.imt_altitude = config.getfloat("SS_MLEO", "imt_altitude")
        self.ss_mleo.imt_lat_deg = config.getfloat("SS_MLEO", "imt_lat_deg")
        self.ss_mleo.imt_long_diff_deg = config.getfloat("SS_MLEO", "imt_long_diff_deg")
        self.ss_mleo.season = config.get("SS_MLEO", "season")
        self.ss_mleo.channel_model = config.get("SS_MLEO", "channel_model")
        self.ss_mleo.antenna_l_s = config.getfloat("SS_MLEO", "antenna_l_s")
        self.ss_mleo.antenna_3_dB = config.getfloat("SS_MLEO", "antenna_3_dB")
        self.ss_mleo.BOLTZMANN_CONSTANT = config.getfloat("SS_MLEO", "BOLTZMANN_CONSTANT")
        self.ss_mleo.EARTH_RADIUS = config.getfloat("SS_MLEO", "EARTH_RADIUS")

        #######################################################################
        # FSS space station
        #######################################################################
        self.fss_ss.frequency               = config.getfloat("FSS_SS", "frequency")
        self.fss_ss.bandwidth               = config.getfloat("FSS_SS", "bandwidth")
        self.fss_ss.tx_power_density        = config.getfloat("FSS_SS", "tx_power_density")
        self.fss_ss.altitude                = config.getfloat("FSS_SS", "altitude")
        self.fss_ss.lat_deg                 = config.getfloat("FSS_SS", "lat_deg")
        self.fss_ss.elevation               = config.getfloat("FSS_SS", "elevation")
        self.fss_ss.azimuth                 = config.getfloat("FSS_SS", "azimuth")
        self.fss_ss.noise_temperature       = config.getfloat("FSS_SS", "noise_temperature")
        self.fss_ss.adjacent_ch_selectivity = config.getfloat("FSS_SS", "adjacent_ch_selectivity")
        self.fss_ss.antenna_gain            = config.getfloat("FSS_SS", "antenna_gain")
        self.fss_ss.antenna_pattern         = config.get("FSS_SS", "antenna_pattern")
        self.fss_ss.imt_altitude            = config.getfloat("FSS_SS", "imt_altitude")
        self.fss_ss.imt_lat_deg             = config.getfloat("FSS_SS", "imt_lat_deg")
        self.fss_ss.imt_long_diff_deg       = config.getfloat("FSS_SS", "imt_long_diff_deg")
        self.fss_ss.season                  = config.get("FSS_SS", "season")
        self.fss_ss.channel_model           = config.get("FSS_SS", "channel_model")
        self.fss_ss.antenna_l_s             = config.getfloat("FSS_SS", "antenna_l_s")
        self.fss_ss.antenna_3_dB            = config.getfloat("FSS_SS", "antenna_3_dB")
        self.fss_ss.BOLTZMANN_CONSTANT      = config.getfloat("FSS_SS", "BOLTZMANN_CONSTANT")
        self.fss_ss.EARTH_RADIUS            = config.getfloat("FSS_SS", "EARTH_RADIUS")

        #######################################################################
        # FSS earth station
        #######################################################################
        self.fss_es.location                = config.get("FSS_ES", "location")
        self.fss_es.x                       = config.getfloat("FSS_ES", "x")
        self.fss_es.y                       = config.getfloat("FSS_ES", "y")
        self.fss_es.min_dist_to_bs          = config.getfloat("FSS_ES", "min_dist_to_bs")
        self.fss_es.max_dist_to_bs          = config.getfloat("FSS_ES", "max_dist_to_bs")
        self.fss_es.height                  = config.getfloat("FSS_ES", "height")
        self.fss_es.elevation_min           = config.getfloat("FSS_ES", "elevation_min")
        self.fss_es.elevation_max           = config.getfloat("FSS_ES", "elevation_max")
        self.fss_es.azimuth                 = config.get("FSS_ES", "azimuth")
        self.fss_es.frequency               = config.getfloat("FSS_ES", "frequency")
        self.fss_es.bandwidth               = config.getfloat("FSS_ES", "bandwidth")
        self.fss_es.adjacent_ch_selectivity = config.getfloat("FSS_ES", "adjacent_ch_selectivity")
        self.fss_es.tx_power_density        = config.getfloat("FSS_ES", "tx_power_density")
        self.fss_es.noise_temperature       = config.getfloat("FSS_ES", "noise_temperature")
        self.fss_es.antenna_gain            = config.getfloat("FSS_ES", "antenna_gain")
        self.fss_es.antenna_pattern         = config.get("FSS_ES", "antenna_pattern")
        self.fss_es.antenna_envelope_gain   = config.getfloat("FSS_ES", "antenna_envelope_gain")
        self.fss_es.diameter                = config.getfloat("FSS_ES", "diameter")
        self.fss_es.channel_model           = config.get("FSS_ES", "channel_model")
        self.fss_es.BOLTZMANN_CONSTANT      = config.getfloat("FSS_ES", "BOLTZMANN_CONSTANT")
        self.fss_es.EARTH_RADIUS            = config.getfloat("FSS_ES", "EARTH_RADIUS")

        # P452 parameters
        self.fss_es.atmospheric_pressure    = config.getfloat("FSS_ES", "atmospheric_pressure")
        self.fss_es.air_temperature         = config.getfloat("FSS_ES", "air_temperature")
        self.fss_es.N0                      = config.getfloat("FSS_ES", "N0")
        self.fss_es.delta_N                 = config.getfloat("FSS_ES", "delta_N")
        self.fss_es.percentage_p            = config.get("FSS_ES", "percentage_p")
        self.fss_es.Dct                     = config.getfloat("FSS_ES", "Dct")
        self.fss_es.Dcr                     = config.getfloat("FSS_ES", "Dcr")
        self.fss_es.Hte                     = config.getfloat("FSS_ES", "Hte")
        self.fss_es.Hre                     = config.getfloat("FSS_ES", "Hre")
        self.fss_es.tx_lat                  = config.getfloat("FSS_ES", "tx_lat")
        self.fss_es.rx_lat                  = config.getfloat("FSS_ES", "rx_lat")
        self.fss_es.polarization            = config.get("FSS_ES", "polarization")
        self.fss_es.clutter_loss            = config.getboolean("FSS_ES", "clutter_loss")

        # HDFSS propagation parameters
        self.fss_es.es_position                  = config.get("FSS_ES", "es_position")
        self.fss_es.shadow_enabled               = config.getboolean("FSS_ES", "shadow_enabled")
        self.fss_es.building_loss_enabled        = config.getboolean("FSS_ES", "building_loss_enabled")
        self.fss_es.same_building_enabled        = config.getboolean("FSS_ES", "same_building_enabled")
        self.fss_es.diffraction_enabled          = config.getboolean("FSS_ES", "diffraction_enabled")
        self.fss_es.bs_building_entry_loss_type  = config.get("FSS_ES", "bs_building_entry_loss_type")
        self.fss_es.bs_building_entry_loss_prob  = config.getfloat("FSS_ES", "bs_building_entry_loss_prob")
        self.fss_es.bs_building_entry_loss_value = config.getfloat("FSS_ES", "bs_building_entry_loss_value")

        #######################################################################
        # Fixed wireless service
        #######################################################################
        self.fs.x                       = config.getfloat("FS", "x")
        self.fs.y                       = config.getfloat("FS", "y")
        self.fs.height                  = config.getfloat("FS", "height")
        self.fs.elevation               = config.getfloat("FS", "elevation")
        self.fs.azimuth                 = config.getfloat("FS", "azimuth")
        self.fs.distribution_enable     = config.get ("FS", "distribution_enable")
        self.fs.distribution_type       = config.get ("FS", "distribution_type")
        self.fs.azimuth_distribution    = config.get ("FS", "azimuth_distribution")
        self.fs.elevation_distribution  = config.get ("FS", "elevation_distribution")
        self.fs.frequency               = config.getfloat("FS", "frequency")
        self.fs.bandwidth               = config.getfloat("FS", "bandwidth")
        self.fs.noise_temperature       = config.getfloat("FS", "noise_temperature")
        self.fs.adjacent_ch_selectivity = config.getfloat("FS", "adjacent_ch_selectivity")
        self.fs.tx_power_density        = config.getfloat("FS", "tx_power_density")
        self.fs.antenna_gain            = config.getfloat("FS", "antenna_gain")
        self.fs.antenna_pattern         = config.get("FS", "antenna_pattern")
        self.fs.diameter                = config.getfloat("FS", "diameter")
        self.fs.channel_model           = config.get("FS", "channel_model")
        self.fs.BOLTZMANN_CONSTANT      = config.getfloat("FS", "BOLTZMANN_CONSTANT")
        self.fs.EARTH_RADIUS            = config.getfloat("FS", "EARTH_RADIUS")
        self.fs.altitude                = config.getfloat("FS", "altitude")

        # P.619 parameters
        self.fs.hibs_lat_deg            = config.getfloat("RAS", "hibs_lat_deg")
        self.fs.imt_altitude            = config.getfloat("FS", "system_altitude")
        self.fs.imt_lat_deg             = config.getfloat("FS", "system_lat_deg")
        self.fs.imt_long_diff_deg       = config.getfloat("FS", "system_long_diff_deg")
        self.fs.season                  = config.get("FS", "season")

        #######################################################################
        # HAPS (airbone) station
        #######################################################################
        self.haps.frequency               = config.getfloat("HAPS", "frequency")
        self.haps.bandwidth               = config.getfloat("HAPS", "bandwidth")
        self.haps.antenna_gain            = config.getfloat("HAPS", "antenna_gain")
        self.haps.tx_power_density        = config.getfloat("HAPS", "eirp_density") - self.haps.antenna_gain - 60
        self.haps.altitude                = config.getfloat("HAPS", "altitude")
        self.haps.lat_deg                 = config.getfloat("HAPS", "lat_deg")
        self.haps.elevation               = config.getfloat("HAPS", "elevation")
        self.haps.azimuth                 = config.getfloat("HAPS", "azimuth")
        self.haps.antenna_pattern         = config.get("HAPS", "antenna_pattern")
        self.haps.imt_altitude            = config.getfloat("HAPS", "imt_altitude")
        self.haps.imt_lat_deg             = config.getfloat("HAPS", "imt_lat_deg")
        self.haps.imt_long_diff_deg       = config.getfloat("HAPS", "imt_long_diff_deg")
        self.haps.season                  = config.get("HAPS", "season")
        self.haps.acs                     = config.getfloat("HAPS", "acs")
        self.haps.channel_model           = config.get("HAPS", "channel_model")
        self.haps.antenna_l_n             = config.getfloat("HAPS", "antenna_l_n")
        self.haps.BOLTZMANN_CONSTANT      = config.getfloat("HAPS", "BOLTZMANN_CONSTANT")
        self.haps.EARTH_RADIUS            = config.getfloat("HAPS", "EARTH_RADIUS")

        #######################################################################
        # RNS
        #######################################################################
        self.rns.x                  = config.getfloat("RNS", "x")
        self.rns.y                  = config.getfloat("RNS", "y")
        self.rns.altitude           = config.getfloat("RNS", "altitude")
        self.rns.frequency          = config.getfloat("RNS", "frequency")
        self.rns.bandwidth          = config.getfloat("RNS", "bandwidth")
        self.rns.noise_temperature  = config.getfloat("RNS", "noise_temperature")
        self.rns.tx_power_density   = config.getfloat("RNS", "tx_power_density")
        self.rns.antenna_gain       = config.getfloat("RNS", "antenna_gain")
        self.rns.antenna_pattern    = config.get("RNS", "antenna_pattern")
        self.rns.season             = config.get("RNS", "season")
        self.rns.imt_altitude       = config.getfloat("RNS", "imt_altitude")
        self.rns.imt_lat_deg        = config.getfloat("RNS", "imt_lat_deg")
        self.rns.channel_model      = config.get("RNS", "channel_model")
        self.rns.acs                = config.getfloat("RNS", "acs")
        self.rns.BOLTZMANN_CONSTANT = config.getfloat("RNS", "BOLTZMANN_CONSTANT")
        self.rns.EARTH_RADIUS       = config.getfloat("RNS", "EARTH_RADIUS")

        #######################################################################
        # ARNS (Air Surveillance and Metereological Radars)
        #######################################################################
        self.arns.x                      = config.getfloat("ARNS", "x")
        self.arns.y                      = config.getfloat("ARNS", "y")
        self.arns.height                 = config.getfloat("ARNS", "height")
        self.arns.elevation              = config.getfloat("ARNS", "elevation")
        self.arns.azimuth                = config.getfloat("ARNS", "azimuth")
        self.arns.distribution_enable    = config.get("ARNS", "distribution_enable")
        self.arns.distribution_type      = config.get("ARNS", "distribution_type")
        self.arns.azimuth_distribution   = config.get("ARNS", "azimuth_distribution")
        self.arns.elevation_distribution = config.get("ARNS", "elevation_distribution")
        self.arns.frequency              = config.getfloat("ARNS", "frequency")
        self.arns.bandwidth              = config.getfloat("ARNS", "bandwidth")
        self.arns.noise_temperature      = config.getfloat("ARNS", "noise_temperature")
        self.arns.tx_power_density       = config.getfloat("ARNS", "tx_power_density")
        self.arns.antenna_gain           = config.getfloat("ARNS", "antenna_gain")
        self.arns.antenna_pattern        = config.get("ARNS", "antenna_pattern")
        self.arns.channel_model          = config.get("ARNS", "channel_model")
        self.arns.acs                    = config.getfloat("ARNS", "acs")
        self.arns.BOLTZMANN_CONSTANT     = config.getfloat("ARNS", "BOLTZMANN_CONSTANT")
        self.arns.EARTH_RADIUS           = config.getfloat("ARNS", "EARTH_RADIUS")
        self.arns.altitude               = config.getfloat("ARNS", "altitude")

        # P.619 parameters
        self.arns.hibs_lat_deg      = config.getfloat("ARNS", "hibs_lat_deg")
        self.arns.imt_altitude      = config.getfloat("ARNS", "system_altitude")
        self.arns.imt_lat_deg       = config.getfloat("ARNS", "system_lat_deg")
        self.arns.imt_long_diff_deg = config.getfloat("ARNS", "system_long_diff_deg")
        self.arns.season            = config.get("ARNS", "season")

        # Air Traffic Control Radar Cossecant Squared antenna parameters
        self.arns.beamwidth_el         = config.getfloat("ARNS", "beamwidth_el")
        self.arns.beamwidth_az         = config.getfloat("ARNS", "beamwidth_az")
        self.arns.csc2_angle           = config.getfloat("ARNS", "maximum_csc2_angle")
        self.arns.highbeam_csc2        = config.getfloat("ARNS", "highbeam_csc2")

        # Phased Array antenna parameters
        self.arns.element_space        = config.getfloat("ARNS", "element_space")

        #######################################################################
        # RAS station
        #######################################################################
        self.ras.x                          = config.getfloat("RAS", "x")
        self.ras.y                          = config.getfloat("RAS", "y")
        self.ras.height                     = config.getfloat("RAS", "height")
        self.ras.elevation                  = config.getfloat("RAS", "elevation")
        self.ras.azimuth                    = config.getfloat("RAS", "azimuth")
        self.ras.distribution_enable        = config.get("RAS", "distribution_enable")
        self.ras.distribution_type          = config.get("RAS", "distribution_type")
        self.ras.azimuth_distribution       = config.get("RAS", "azimuth_distribution")
        self.ras.elevation_distribution     = config.get("RAS", "elevation_distribution")
        self.ras.frequency                  = config.getfloat("RAS", "frequency")
        self.ras.bandwidth                  = config.getfloat("RAS", "bandwidth")
        self.ras.antenna_noise_temperature  = config.getfloat("RAS", "antenna_noise_temperature")
        self.ras.receiver_noise_temperature = config.getfloat("RAS", "receiver_noise_temperature")
        self.ras.adjacent_ch_selectivity    = config.getfloat("FSS_ES", "adjacent_ch_selectivity")
        self.ras.antenna_efficiency         = config.getfloat("RAS", "antenna_efficiency")
        self.ras.antenna_gain               = config.getfloat("RAS", "antenna_gain")
        self.ras.antenna_pattern            = config.get("RAS", "antenna_pattern")
        self.ras.diameter                   = config.getfloat("RAS", "diameter")
        self.ras.channel_model              = config.get("RAS", "channel_model")
        self.ras.BOLTZMANN_CONSTANT         = config.getfloat("RAS", "BOLTZMANN_CONSTANT")
        self.ras.EARTH_RADIUS               = config.getfloat("RAS", "EARTH_RADIUS")
        self.ras.SPEED_OF_LIGHT             = config.getfloat("RAS", "SPEED_OF_LIGHT")
        self.ras.altitude                   = config.getfloat("RAS", "altitude")

        # P.619 parameters
        self.ras.hibs_lat_deg               = config.getfloat("RAS", "hibs_lat_deg")
        self.ras.imt_altitude               = config.getfloat("RAS", "system_altitude")
        self.ras.imt_lat_deg                = config.getfloat("RAS", "system_lat_deg")
        self.ras.imt_long_diff_deg          = config.getfloat("RAS", "system_long_diff_deg")
        self.ras.season                     = config.get("RAS", "season")

        # P452 parameters
        self.ras.atmospheric_pressure       = config.getfloat("RAS", "atmospheric_pressure")
        self.ras.air_temperature            = config.getfloat("RAS", "air_temperature")
        self.ras.N0                         = config.getfloat("RAS", "N0")
        self.ras.delta_N                    = config.getfloat("RAS", "delta_N")
        self.ras.percentage_p               = config.get("RAS", "percentage_p")
        self.ras.Dct                        = config.getfloat("RAS", "Dct")
        self.ras.Dcr                        = config.getfloat("RAS", "Dcr")
        self.ras.Hte                        = config.getfloat("RAS", "Hte")
        self.ras.Hre                        = config.getfloat("RAS", "Hre")
        self.ras.tx_lat                     = config.getfloat("RAS", "tx_lat")
        self.ras.rx_lat                     = config.getfloat("RAS", "rx_lat")
        self.ras.polarization               = config.get("RAS", "polarization")
        self.ras.clutter_loss               = config.getboolean("RAS", "clutter_loss")

        #######################################################################
        # EESS passive
        #######################################################################
        self.eess_passive.frequency               = config.getfloat("EESS_PASSIVE", "frequency")
        self.eess_passive.bandwidth               = config.getfloat("EESS_PASSIVE", "bandwidth")
        self.eess_passive.nadir_angle             = config.getfloat("EESS_PASSIVE", "nadir_angle")
        self.eess_passive.altitude                = config.getfloat("EESS_PASSIVE", "altitude")
        self.eess_passive.antenna_pattern         = config.get("EESS_PASSIVE", "antenna_pattern")
        self.eess_passive.antenna_efficiency      = config.getfloat("EESS_PASSIVE", "antenna_efficiency")
        self.eess_passive.antenna_diameter        = config.getfloat("EESS_PASSIVE", "antenna_diameter")
        self.eess_passive.antenna_gain            = config.getfloat("EESS_PASSIVE", "antenna_gain")
        self.eess_passive.channel_model           = config.get("EESS_PASSIVE", "channel_model")
        self.eess_passive.imt_altitude            = config.getfloat("EESS_PASSIVE", "imt_altitude")
        self.eess_passive.imt_lat_deg             = config.getfloat("EESS_PASSIVE", "imt_lat_deg")
        self.eess_passive.season                  = config.get("EESS_PASSIVE", "season")
        self.eess_passive.BOLTZMANN_CONSTANT      = config.getfloat("EESS_PASSIVE", "BOLTZMANN_CONSTANT")
        self.eess_passive.EARTH_RADIUS            = config.getfloat("EESS_PASSIVE", "EARTH_RADIUS")

        #######################################################################
        # IMT BASE STATION
        #######################################################################
        self.imt_bs.topology                      = config.get("IMT BASE STATION", "topology")
        self.imt_bs.wrap_around                   = config.getboolean("IMT BASE STATION", "wrap_around")
        self.imt_bs.num_clusters                  = config.getint("IMT BASE STATION", "num_clusters")
        self.imt_bs.intersite_distance            = config.getfloat("IMT BASE STATION", "intersite_distance")
        self.imt_bs.frequency                     = config.getfloat("IMT BASE STATION", "frequency")
        self.imt_bs.bandwidth                     = config.getfloat("IMT BASE STATION", "bandwidth")
        self.imt_bs.rb_bandwidth                  = config.getfloat("IMT BASE STATION", "rb_bandwidth")
        self.imt_bs.spectral_mask                 = config.get("IMT BASE STATION", "spectral_mask")
        self.imt_bs.spurious_emissions            = config.getfloat("IMT BASE STATION", "spurious_emissions")
        self.imt_bs.guard_band_ratio              = config.getfloat("IMT BASE STATION", "guard_band_ratio")
        self.imt_bs.bs_load_probability           = config.getfloat("IMT BASE STATION", "bs_load_probability")
        self.imt_bs.bs_height                     = config.getfloat("IMT BASE STATION", "bs_height")
        self.imt_bs.bs_noise_figure               = config.getfloat("IMT BASE STATION", "bs_noise_figure")
        self.imt_bs.bs_noise_temperature          = config.getfloat("IMT BASE STATION", "bs_noise_temperature")
        self.imt_bs.bs_ohmic_loss                 = config.getfloat("IMT BASE STATION", "bs_ohmic_loss")
        self.imt_bs.dl_attenuation_factor         = config.getfloat("IMT BASE STATION", "dl_attenuation_factor")
        self.imt_bs.dl_sinr_min                   = config.getfloat("IMT BASE STATION", "dl_sinr_min")
        self.imt_bs.dl_sinr_max                   = config.getfloat("IMT BASE STATION", "dl_sinr_max")
        self.imt_bs.channel_model                 = config.get("IMT BASE STATION", "channel_model")
        self.imt_bs.los_adjustment_factor         = config.getfloat("IMT BASE STATION", "los_adjustment_factor")
        self.imt_bs.shadowing                     = config.getboolean("IMT BASE STATION", "shadowing")
        self.imt_bs.noise_temperature             = config.getfloat("IMT BASE STATION", "noise_temperature")
        self.imt_bs.BOLTZMANN_CONSTANT            = config.getfloat("IMT BASE STATION", "BOLTZMANN_CONSTANT")
        self.imt_bs.EARTH_RADIUS                  = config.getfloat("IMT BASE STATION", "EARTH_RADIUS")
        self.imt_bs.altitude                      = config.getfloat("IMT BASE STATION", "altitude")

        # P.619 parameters
        self.imt_bs.hibs_lat_deg            = config.getfloat("IMT BASE STATION", "hibs_lat_deg")
        self.imt_bs.imt_altitude            = config.getfloat("IMT BASE STATION", "system_altitude")
        self.imt_bs.imt_lat_deg             = config.getfloat("IMT BASE STATION", "system_lat_deg")
        self.imt_bs.imt_long_diff_deg       = config.getfloat("IMT BASE STATION", "system_long_diff_deg")
        self.imt_bs.season                  = config.get("IMT BASE STATION", "season")

        # IMT Base Station Antenna parameters

        self.antenna_imt_bs.bs_antenna_type            = config.get("IMT BASE STATION", "bs_antenna_type")
        self.antenna_imt_bs.adjacent_antenna_model     = config.get("IMT BASE STATION", "adjacent_antenna_model")
        self.antenna_imt_bs.bs_normalization           = config.getboolean("IMT BASE STATION", "bs_normalization")
        self.antenna_imt_bs.bs_normalization_file      = config.get("IMT BASE STATION", "bs_normalization_file")
        self.antenna_imt_bs.bs_element_pattern         = config.get("IMT BASE STATION", "bs_element_pattern")

        self.antenna_imt_bs.bs_element_max_g           = config.getfloat("IMT BASE STATION", "bs_element_max_g")
        self.antenna_imt_bs.bs_element_phi_3db         = config.getfloat("IMT BASE STATION", "bs_element_phi_3db")
        self.antenna_imt_bs.bs_element_theta_3db       = config.getfloat("IMT BASE STATION", "bs_element_theta_3db")
        self.antenna_imt_bs.bs_element_am              = config.getfloat("IMT BASE STATION", "bs_element_am")
        self.antenna_imt_bs.bs_element_sla_v           = config.getfloat("IMT BASE STATION", "bs_element_sla_v")
        self.antenna_imt_bs.bs_n_rows                  = config.getfloat("IMT BASE STATION", "bs_n_rows")
        self.antenna_imt_bs.bs_n_columns               = config.getfloat("IMT BASE STATION", "bs_n_columns")
        self.antenna_imt_bs.bs_element_horiz_spacing   = config.getfloat("IMT BASE STATION", "bs_element_horiz_spacing")
        self.antenna_imt_bs.bs_element_vert_spacing    = config.getfloat("IMT BASE STATION", "bs_element_vert_spacing")
        self.antenna_imt_bs.bs_multiplication_factor   = config.getfloat("IMT BASE STATION", "bs_multiplication_factor")
        self.antenna_imt_bs.bs_minimum_array_gain      = config.getfloat("IMT BASE STATION", "bs_minimum_array_gain")
        self.antenna_imt_bs.bs_downtilt                = config.getfloat("IMT BASE STATION", "bs_downtilt")
        self.antenna_imt_bs.bf_enable                      = config.get("IMT BASE STATION", "bf_enable")
