import numpy as np
from dataclasses import dataclass, field, asdict
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_orbit import ParametersOrbit
from sharc.parameters.imt.parameters_imt_mss_dc import ParametersSelectActiveSatellite, ParametersSectorPositioning
from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528


@dataclass
class ParametersMssD2d(ParametersBase):
    """Define parameters for a MSS-D2D - NGSO Constellation."""
    section_name: str = "mss_d2d"

    nested_parameters_enabled: bool = True

    is_space_to_earth: bool = True

    # MSS_D2D system name
    name: str = "Default"

    # Orbit parameters
    orbits: list[ParametersOrbit] = field(default_factory=lambda: [ParametersOrbit()])

    # MSS_D2D system center frequency in MHz
    frequency: float = 2110.0

    # MSS_D2d system bandwidth in MHz
    bandwidth: float = 5.0

    # In case you want to use a load factor for beams
    # that means that each beam has a probability of `beams_load_factor` to be active
    beams_load_factor: float = 1.0

    # Central beam positioning
    center_beam_positioning: ParametersSectorPositioning = field(default_factory=ParametersSectorPositioning)

    # Adjacent channel emissions type
    # Possible values are "ACLR", "SPECTRAL_MASK" and "OFF"
    adjacent_ch_emissions: str = "OFF"

    # Transmitter spectral mask
    spectral_mask: str = "MSS"

    # Out-of-band spurious emissions in dB/MHz
    spurious_emissions: float = -13.0

    # Adjacent channel leakage ratio in dB
    adjacent_ch_leak_ratio: float = 45.0

    # Single beam cell radius in network topology [m]
    cell_radius: float = 19000.0

    # NTN Intersite Distance (m). Intersite distance = Cell Radius * sqrt(3)
    intersite_distance: float = np.sqrt(3) * cell_radius

    # Satellite tx power density in dBW/MHz
    tx_power_density: float = 40.0

    # # Satellite Tx max Gain in dBi
    # antenna_gain: float = 30.0

    # Number of beams per satellite
    num_sectors: int = 19

    # Satellite antenna pattern
    # Antenna pattern from ITU-R S.1528
    # Possible values: "ITU-R-S.1528-Taylor", "ITU-R-S.1528-LEO"
    antenna_pattern: str = "ITU-R-S.1528-Taylor"

    # Radius of the antenna's circular aperture in meters
    antenna_diamter: float = 1.0

    # The required near-in-side-lobe level (dB) relative to peak gain
    antenna_l_s: float = -6.75

    # 3 dB beamwidth angle (3 dB below maximum gain) [degrees]
    antenna_3_dB_bw: float = 4.4127

    # Paramters for the ITU-R-S.1528 antenna patterns
    antenna_s1528: ParametersAntennaS1528 = field(default_factory=ParametersAntennaS1528)

    sat_is_active_if: ParametersSelectActiveSatellite = field(default_factory=ParametersSelectActiveSatellite)

    # paramters for channel model
    param_p619: ParametersP619 = field(default_factory=ParametersP619)
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

        self.propagate_parameters()

        self.validate(self.section_name)

    def __post_init__(self):
        self.beam_radius = self.cell_radius
        self.num_beams = self.num_sectors

    def validate(self, ctx):
        # Now do the sanity check for some parameters
        if self.num_sectors not in [1, 7, 19]:
            raise ValueError(f"ParametersMssD2d: Invalid number of sectors {self.num_sectors}")
        self.num_beams = self.num_sectors

        if self.cell_radius <= 0:
            raise ValueError(f"ParametersMssD2d: cell_radius must be greater than 0, but is {self.cell_radius}")
        else:
            self.intersite_distance = np.sqrt(3) * self.cell_radius
            self.beam_radius = self.cell_radius

        if self.adjacent_ch_emissions not in ["SPECTRAL_MASK", "ACLR", "OFF"]:
            raise ValueError(f"""ParametersMssD2d: Invalid adjacent channel emissions {self.adjacent_ch_emissions}""")

        if self.spectral_mask.upper() not in ["IMT-2020", "3GPP E-UTRA", "MSS"]:
            raise ValueError(f"""ParametersMssD2d: Inavlid Spectral Mask Name {self.spectral_mask}""")

        if self.channel_model.upper() not in ["FSPL", "P619", "SATELLITESIMPLE"]:
            raise ValueError(f"Invalid channel model name {self.channel_model}")

        if self.beams_load_factor < 0.0 or self.beams_load_factor > 1.0:
            raise ValueError(f"{ctx}.beams_load_factor must be in interval [0.0, 1.0]")

        super().validate(ctx)

    def propagate_parameters(self):
        self.antenna_s1528.set_external_parameters(antenna_pattern=self.antenna_pattern,
                                                   frequency=self.frequency,
                                                   bandwidth=self.bandwidth,
                                                   antenna_l_s=self.antenna_l_s,
                                                   antenna_3_dB_bw=self.antenna_3_dB_bw,)

        if self.channel_model == "P619":
            # mean station altitude in meters
            m_alt = 0
            for orbit in self.orbits:
                m_alt += orbit.perigee_alt_km * 1e3
            m_alt /= len(self.orbits)

            self.param_p619.set_external_parameters(
                space_station_alt_m=m_alt,
                earth_station_alt_m=self.earth_station_alt_m,
                earth_station_lat_deg=self.earth_station_lat_deg,
                earth_station_long_diff_deg=self.earth_station_long_diff_deg,
                season=self.season
            )


if __name__ == "__main__":
    # Run validation for input parameters
    import os
    import pprint

    # Load default simulator parameters
    yaml_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../input/parameters.yaml")
    mss_d2d_params = ParametersMssD2d()
    mss_d2d_params.load_parameters_from_file(yaml_file_path)
    pprint.pprint(asdict(mss_d2d_params), sort_dicts=False)
