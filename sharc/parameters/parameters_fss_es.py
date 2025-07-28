# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sharc.support.sharc_utils import is_float
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_p452 import ParametersP452
from sharc.parameters.parameters_p619 import ParametersP619


@dataclass
class ParametersFssEs(ParametersBase):
    """Dataclass containing the Fixed Satellite Services - Earth Station
    parameters for the simulator
    """
    section_name: str = "fss_es"

    # type of FSS-ES location:
    # FIXED - position must be given
    # CELL - random within central cell
    # NETWORK - random within whole network
    # UNIFORM_DIST - uniform distance from cluster centre,
    #                between min_dist_to_bs and max_dist_to_bs
    location: str = "UNIFORM_DIST"
    # x-y coordinates [m] (only if FIXED location is chosen)
    x: float = 10000.0
    y: float = 0.0
    # minimum distance from BSs [m]
    min_dist_to_bs: float = 10.0
    # maximum distance from centre BSs [m] (only if UNIFORM_DIST is chosen)
    max_dist_to_bs: float = 10.0
    # antenna height [m]
    height: float = 6.0
    # Elevation angle [deg], minimum and maximum, values
    elevation_min: float = 48.0
    elevation_max: float = 80.0
    # Azimuth angle [deg]
    # either a specific angle or string 'RANDOM'
    azimuth: str = 0.2
    # center frequency [MHz]
    frequency: float = 43000.0
    # bandwidth [MHz]
    bandwidth: float = 6.0
    # adjacent channel selectivity (dB)
    adjacent_ch_selectivity: float = 0.0
    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density: float = -68.3
    # System receive noise temperature [K]
    noise_temperature: float = 950.0
    # antenna peak gain [dBi]
    antenna_gain: float = 0.0
    # Antenna pattern of the FSS Earth station
    # Possible values: "ITU-R S.1855", "ITU-R S.465", "ITU-R S.580", "OMNI",
    #                  "Modified ITU-R S.465"
    antenna_pattern: str = "Modified ITU-R S.465"
    # Antenna envelope gain (dBi) - only relevant for "Modified ITU-R S.465" model
    antenna_envelope_gain: float = 0.0
    # Diameter of the antenna [m]
    diameter: float = 1.8
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "TerrestrialSimple" (FSPL + clutter loss)
    #                                    "P452"
    #                                    "TVRO-URBAN"
    #                                    "TVRO-SUBURBAN"
    #                                    "HDFSS"
    channel_model: str = "P452"

    # P452 parameters
    param_p452 = ParametersP452()
    # Total air pressure in hPa
    atmospheric_pressure: float = 935.0
    # Temperature in Kelvin
    air_temperature: float = 300.0
    # Sea-level surface refractivity (use the map)
    N0: float = 352.58
    # Average radio-refractive (use the map)
    delta_N: float = 43.127
    # Percentage p. Float (0 to 100) or RANDOM
    percentage_p: str = "0.2"
    # Distance over land from the transmit and receive antennas to the coast (km)
    Dct: float = 70.0
    # Distance over land from the transmit and receive antennas to the coast (km)
    Dcr: float = 70.0
    # Effective height of interfering antenna (m)
    Hte: float = 20.0
    # Effective height of interfered-with antenna (m)
    Hre: float = 3.0
    # Latitude of transmitter
    tx_lat: float = -23.55028
    # Latitude of receiver
    rx_lat: float = -23.17889
    # Antenna polarization
    polarization: str = "horizontal"
    # Determine whether clutter loss following ITU-R P.2108 is added (TRUE/FALSE)
    clutter_loss: bool = True
    clutter_type: str = "one-end"

    # Parameters for the P.619 propagation model used for sharing studies between IMT-NTN and FSS-ES
    #    space_station_alt_m - altiteude of the IMT-MSS station
    #    earth_station_alt_m - altitude of FSS-ES system (in meters)
    #    earth_station_lat_deg - latitude of FSS-ES system (in degrees)
    #    earth_station_long_diff_deg - difference between longitudes of IMT-NTN station and FSS-ES system
    #      (positive if space-station is to the East of earth-station)
    #    season - season of the year.
    param_p619 = ParametersP619()
    space_station_alt_m: float = 35780000.0
    earth_station_alt_m: float = 0.0
    earth_station_lat_deg: float = 0.0
    earth_station_long_diff_deg: float = 0.0
    season: str = "SUMMER"

    # HDFSS propagation parameters
    # HDFSS position relative to building it is on. Possible values are
    # ROOFTOP and BUILDINGSIDE
    es_position: str = "ROOFTOP"
    # Enable shadowing loss
    shadow_enabled: bool = True
    # Enable building entry loss
    building_loss_enabled: bool = True
    # Enable interference from IMT stations at the same building as the HDFSS
    same_building_enabled: bool = False
    # Enable diffraction loss
    diffraction_enabled: bool = False
    # Building entry loss type applied between BSs and HDFSS ES. Options are:
    # P2109_RANDOM: random probability at P.2109 model, considering elevation
    # P2109_FIXED: fixed probability at P.2109 model, considering elevation.
    #              Probability must be specified in bs_building_entry_loss_prob.
    # FIXED_VALUE: fixed value per BS. Value must be specified in
    #              bs_building_entry_loss_value.
    bs_building_entry_loss_type: str = "P2109_FIXED"
    # Probability of building entry loss not exceeded if
    # bs_building_entry_loss_type = P2109_FIXED
    bs_building_entry_loss_prob: float = 0.75
    # Value in dB of building entry loss if
    # bs_building_entry_loss_type = FIXED_VALUE
    bs_building_entry_loss_value: float = 35

    def load_parameters_from_file(self, config_file: str):
        """Load the parameters from file an run a sanity check

        Parameters
        ----------
        file_name : str
            the path to the configuration file

        Raises
        ------
        ValueError
            if a parameter is not valid
        """
        super().load_parameters_from_file(config_file)

        if self.location not in ["FIXED", "CELL", "NETWORK", "UNIFORM_DIST"]:
            raise ValueError(f"ParametersFssEs: \
                             Invalid value for paramter location - {self.location}. \
                            Allowed values are \"FIXED\", \"CELL\", \"NETWORK\", \"UNIFORM_DIST\".")

        if self.antenna_pattern not in [
            "ITU-R S.1855", "ITU-R S.465", "ITU-R S.580", "OMNI",
            "Modified ITU-R S.465",
        ]:
            raise ValueError(f"ParametersFssEs: \
                             Invalid value for paramter antenna_pattern - {self.antenna_pattern}. \
                             Allowed values are \
                             \"ITU-R S.1855\", \"ITU-R S.465\", \"ITU-R S.580\", \"OMNI\", \
                             \"Modified ITU-R S.465\"")

        if is_float(self.azimuth):
            self.azimuth = float(self.azimuth)
        elif self.azimuth.upper() != "RANDOM":
            if self.azimuth.isnumeric():
                self.azimuth = float(self.azimuth)
            else:
                raise ValueError(f"""ParametersFssEs:
                                Invalid value for parameter azimuth - {self.azimuth}.
                                Allowed values are \"RANDOM\" or a angle in degrees.""")

        if isinstance(self.percentage_p, str) and self.percentage_p.upper() != "RANDOM":
            percentage_p = float(self.percentage_p)
            if percentage_p <= 0 or percentage_p > 1:
                raise ValueError(f"""ParametersFssEs:
                        Invalid value for parameter percentage_p - {self.percentage_p}.
                        Allowed values are a percentage between 0 and 1 (exclusive).""")

        if self.polarization.lower() not in ["horizontal", "vertical"]:
            raise ValueError(f"ParametersFssEss: \
                             Invalid value for parameter polarization - {self.polarization}. \
                             Allowed values are: \"horizontal\", \"vertical\"")

        if self.es_position.upper() not in ["BUILDINGSIDE", "ROOFTOP"]:
            raise ValueError(f"ParametersFssEss: \
                             Invalid value for parameter es_position - {self.es_position} \
                             Allowed values are \"BUILDINGSIDE\", \"ROOFTOP\".")

        if self.bs_building_entry_loss_type not in ["P2109_RANDOM", "P2109_FIXED", "FIXED_VALUE"]:
            raise ValueError(f"ParametersFssEs: \
                             Invalid value for parameter bs_building_entry_loss_type - \
                             {self.bs_building_entry_loss_type} \
                             Allowd values are \"P2109_RANDOM\", \"P2109_FIXED\", \"FIXED_VALUE\".")
        if self.channel_model.upper() not in [
            "FSPL", "TERRESTRIALSIMPLE", "P452", "P619",
            "TVRO-URBAN", "TVRO-SUBURBAN", "HDFSS", "UMA", "UMI",
        ]:
            raise ValueError(
                f"ParametersFssEs: Invalid value for parameter channel_model - {self.channel_model}",
            )

        if self.channel_model == "P452":
            self.param_p452.load_from_paramters(self)

        elif self.channel_model == "P619":
            self.param_p619.load_from_paramters(self)
