from dataclasses import dataclass
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_p452 import ParametersP452


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
    # "TerrestrialSimple" (FSPL + clutter loss),
    # P452
    channel_model: str = "FSPL"

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
        if self.channel_model not in ["FSPL", "TerrestrialSimple", "P452"]:
            raise ValueError(
                "Invalid channel_model, must be either 'FSPL', 'TerrestrialSimple', or 'P452'",
            )
        if self.channel_model == "P452":
            self.param_p452.load_from_paramters(self)
