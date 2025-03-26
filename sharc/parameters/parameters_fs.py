from dataclasses import dataclass
from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersFs(ParametersBase):
    """
    Parameters definitions for fixed wireless service systems.
    """

    section_name: str = "fs"

    # x-y coordinates [meters]
    x: float = 1000.0
    y: float = 0.0

    # Antenna height [meters]
    height: float = 15.0

    # Elevation angle [degrees]
    elevation: float = -10.0

    # Azimuth angle [degrees]
    azimuth: float = 180.0

    # Center frequency [MHz]
    frequency: float = 27250.0

    # Bandwidth [MHz]
    bandwidth: float = 112.0

    # System receive noise temperature [Kelvin]
    noise_temperature: float = 290.0

    # Adjacent channel selectivity [dB]
    adjacent_ch_selectivity: float = 20.0

    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density: float = -68.3

    # Antenna peak gain [dBi]
    antenna_gain: float = 36.9

    # Antenna pattern of the fixed wireless service
    # Possible values: "ITU-R F.699", "OMNI"
    antenna_pattern: str = "ITU-R F.699"

    # Diameter of antenna [meters]
    diameter: float = 0.3

    # Channel model, possible values are "FSPL" (free-space path loss),
    # "TerrestrialSimple" (FSPL + clutter loss)
    channel_model: str = "FSPL"

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

        # Implementing sanity checks for critical parameters
        if not (-90 <= self.elevation <= 90):
            raise ValueError(
                "Elevation angle must be between -90 and 90 degrees.",
            )

        if not (0 <= self.azimuth <= 360):
            raise ValueError(
                "Azimuth angle must be between 0 and 360 degrees.",
            )

        if self.antenna_pattern not in ["ITU-R F.699", "OMNI"]:
            raise ValueError(
                f"Invalid antenna_pattern: {self.antenna_pattern}",
            )

        # Sanity check for channel model
        if self.channel_model not in ["FSPL", "TerrestrialSimple"]:
            raise ValueError(
                "Invalid channel_model, must be either 'FSPL' or 'TerrestrialSimple'",
            )
