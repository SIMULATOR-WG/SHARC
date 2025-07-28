from dataclasses import dataclass

from sharc.parameters.parameters_space_station import ParametersSpaceStation

# The default values come from Report ITU-R SA.2488-0, table 19 (the only earth-to-space MetSat entry)
# TODO: let MetSat as interferrer
# TODO: ver com professor se considerar as tabelas do report está correto


@dataclass
class ParametersMetSatSS(ParametersSpaceStation):
    """
    Defines parameters for MetSat space stations (SS)
    and their interaction with other services based on ITU recommendations.
    """
    section_name: str = "mestat_ss"

    # raw data transmission in 8175-8215 range
    frequency: float = 8195.0  # Satellite center frequency [MHz]
    bandwidth: float = 20.0

    # Elevation angle [deg]
    elevation: float = 3  # Minimum == 3

    # Satellite altitude [m]
    # Altitude of the satellite above the Earth's surface in meters
    altitude: float = 35786000.0  # 35786000 for GSO

    # Antenna pattern of the satellite
    # Possible values: "ITU-R S.672"
    antenna_pattern: str = "ITU-R S.672"  # Antenna radiation pattern

    # antenna peak gain [dBi]
    antenna_gain: float = 52.0

    # The required near-in-side-lobe level (dB) relative to peak gain
    # according to ITU-R S.672-4
    # TODO: check if this changes from fss_ss
    antenna_l_s: float = -20.0
    # 3 dB beamwidth angle (3 dB below maximum gain) [degrees]
    antenna_3_dB: float = 0.65

    # Channel model, possible values are "FSPL" (free-space path loss), "P619"
    # See params space station to check P619 needed params
    channel_model: str = "FSPL"  # Channel model to be used

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
        print(self.antenna_pattern)
        if self.antenna_pattern not in ["ITU-R S.672"]:
            raise ValueError(f"Invalid antenna_pattern: {self.antenna_pattern}")

        # Check channel model
        if self.channel_model not in ["FSPL", "P619"]:
            raise ValueError(
                "Invalid channel_model, must be either 'FSPL' or 'P619'",
            )

        # Check season
        if self.season not in ["SUMMER", "WINTER"]:
            raise ValueError(
                "Invalid season, must be either 'SUMMER' or 'WINTER'",
            )
