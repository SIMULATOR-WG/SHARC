# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersRns(ParametersBase):
    """
    Simulation parameters for radionavigation service
    """
    section_name: str = "rns"
    # x-y coordinates [m]
    x: float = 660.0
    y: float = -370.0
    # altitude [m]
    altitude: float = 150.0
    # center frequency [MHz]
    frequency: float = 32000.0
    # bandwidth [MHz]
    bandwidth: float = 60.0
    # System receive noise temperature [K]
    noise_temperature: float = 1154.0
    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density: float = -70.79
    # antenna peak gain [dBi]
    antenna_gain: float = 30.0
    # Antenna pattern of the fixed wireless service
    # Possible values: "ITU-R M.1466", "OMNI"
    antenna_pattern: str = "ITU-R M.1466"
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "SatelliteSimple" (FSPL + 4 dB + clutter loss)
    #                                    "P619"
    channel_model: str = "P619"
    # Parameters for the P.619 propagation model
    #    earth_station_alt_m - altitude of IMT system (in meters)
    #    earth_station_lat_deg - latitude of IMT system (in degrees)
    #    earth_station_long_diff_deg - difference between longitudes of IMT and satellite system
    #      (positive if space-station is to the East of earth-station)
    #    season - season of the year.
    earth_station_alt_m: float = 0.0
    earth_station_lat_deg: float = 0.0
    earth_station_long_diff_deg: float = 0.0
    season: str = "SUMMER"
    # Adjacent channel selectivity [dB]
    acs: float = 30.0

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
        if self.antenna_pattern not in ["ITU-R M.1466", "OMNI"]:
            raise ValueError(f"ParametersRns: \
                             Invalid value for parameter {self.antenna_pattern}. \
                             Allowed values are \"ITU-R M.1466\", \"OMNI\".")
        if self.channel_model.upper() not in ["FSPL", "SatelliteSimple", "P619"]:
            raise ValueError(f"ParametersRns: \
                             Invalid value for paramter channel_model = {self.channel_model}. \
                             Possible values are \"FSPL\", \"SatelliteSimple\", \"P619\".")
        if self.season.upper() not in ["SUMMER", "WINTER"]:
            raise ValueError(f"ParametersRns: \
                             Invalid value for parameter season - {self.season}. \
                             Possible values are \"SUMMER\", \"WINTER\".")
