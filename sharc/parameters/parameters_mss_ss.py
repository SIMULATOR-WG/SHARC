# -*- coding: utf-8 -*-
"""Parameters definitions for NTN systems
"""
from dataclasses import dataclass, field

from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528
import numpy as np


@dataclass
class ParametersMssSs(ParametersBase):
    """
    Simulation parameters for Mobile Satellite System - Space Station.
    """
    section_name: str = "mss_ss"

    is_space_to_earth: bool = True

    # Satellite bore-sight - this is the center of the central beam.
    x: float = 0.0
    y: float = 0.0

    # MSS_SS system center frequency in MHz
    frequency: float = 2110.0

    # MSS_SS system bandwidth in MHz
    bandwidth: float = 20.0

    # Transmitter spectral mask
    spectral_mask: str = "3GPP E-UTRA"

    # Out-of-band spurious emissions in dB/MHz
    spurious_emissions: float = -13

    # Adjacent channel leakage ratio in dB
    adjacent_ch_leak_ratio: float = 45.0

    # MSS_SS altitude w.r.t. sea level
    altitude: float = 1200000.0

    # NTN cell radius in network topology [m]
    cell_radius: float = 45000.0

    # NTN Intersite Distance (m). Intersite distance = Cell Radius * sqrt(3)
    intersite_distance: float = np.sqrt(3) * cell_radius

    # Satellite tx power density in dBW/MHz
    tx_power_density: float = 40.0

    # Satellite Tx max Gain in dBi
    antenna_gain: float = 30.0

    # Satellite azimuth w.r.t. simulation x-axis
    azimuth: float = 45.0

    # Satellite elevation w.r.t. simulation xy-plane (horizon)
    elevation: float = 90.0

    # Number of sectors
    num_sectors: int = 7

    # Satellite antenna pattern
    # Antenna pattern from ITU-R S.1528
    # Possible values: "ITU-R-S.1528-Section1.2", "ITU-R-S.1528-LEO"
    antenna_pattern: str = "ITU-R-S.1528-LEO"

    # Radius of the antenna's circular aperture in meters
    antenna_diamter: float = 1.0

    # The required near-in-side-lobe level (dB) relative to peak gain
    antenna_l_s: float = -6.75

    # 3 dB beamwidth angle (3 dB below maximum gain) [degrees]
    antenna_3_dB_bw: float = 4.4127

    # Paramters for the ITU-R-S.1528 antenna patterns
    antenna_s1528: ParametersAntennaS1528 = field(default_factory=ParametersAntennaS1528)

    # paramters for channel model
    param_p619: ParametersP619 = field(default_factory=ParametersP619)
    space_station_alt_m: float = 35780000.0
    earth_station_alt_m: float = 0.0
    earth_station_lat_deg: float = 0.0
    earth_station_long_diff_deg: float = 0.0
    season: str = "SUMMER"
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "SatelliteSimple" (FSPL + 4 + clutter loss)
    #                                    "P619"
    channel_model: str = "P619"

    def load_parameters_from_file(self, config_file: str):
        """
        Load the parameters from a file and run a sanity check.

        Parameters
        ----------
        config_file : str
            The path to the configuration file.

        Raises
        ------
        ValueError
            If a parameter is not valid.
        """
        super().load_parameters_from_file(config_file)

        # Now do the sanity check for some parameters
        if self.num_sectors not in [1, 7, 19]:
            raise ValueError(f"ParametersMssSs: Invalid number of sectors {self.num_sectors}")

        if self.cell_radius <= 0:
            raise ValueError(f"ParametersMssSs: cell_radius must be greater than 0, but is {self.cell_radius}")
        else:
            self.intersite_distance = np.sqrt(3) * self.cell_radius

        if not np.all((0 <= self.azimuth) & (self.azimuth <= 360)):
            raise ValueError(
                "ParametersMssSs: bs_azimuth values must be between 0 and 360 degrees")

        if not np.all((0 <= self.elevation) & (self.elevation <= 90)):
            raise ValueError(
                "ParametersMssSs: bs_elevation values must be between 0 and 90 degrees")

        if self.spectral_mask.upper() not in ["IMT-2020", "3GPP E-UTRA"]:
            raise ValueError(f"""ParametersImt: Inavlid Spectral Mask Name {self.spectral_mask}""")

        self.antenna_s1528.set_external_parameters(frequency=self.frequency,
                                                   bandwidth=self.bandwidth,
                                                   antenna_gain=self.antenna_gain,
                                                   antenna_l_s=self.antenna_l_s,
                                                   antenna_3_dB_bw=self.antenna_3_dB_bw,)

        if self.channel_model.upper() not in ["FSPL", "P619", "SATELLITESIMPLE"]:
            raise ValueError(f"Invalid channel model name {self.channel_model}")

        if self.channel_model == "P619":
            self.param_p619.set_external_parameters(
                space_station_alt_m=self.altitude,
                earth_station_alt_m=self.earth_station_alt_m,
                earth_station_lat_deg=self.earth_station_lat_deg,
                earth_station_long_diff_deg=self.earth_station_long_diff_deg,
                season=self.season
            )
