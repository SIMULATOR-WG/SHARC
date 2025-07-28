from dataclasses import dataclass

from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.parameters_space_station import ParametersSpaceStation


@dataclass
class ParametersEessSS(ParametersSpaceStation):
    """
    Defines parameters for Earth Exploration Satellite Service (EESS) sensors/space stations (SS)
    and their interaction with other services based on ITU recommendations.
    """
    section_name: str = "eess_ss"

    is_space_to_earth: bool = True

    # Sensor center frequency [MHz]
    frequency: float = 23900.0  # Center frequency of the sensor in MHz

    # Sensor bandwidth [MHz]
    bandwidth: float = 200.0  # Bandwidth of the sensor in MHz

    # Off-nadir pointing angle [deg]
    nadir_angle: float = 46.6  # Angle in degrees away from the nadir point

    # Sensor altitude [m]
    # Altitude of the sensor above the Earth's surface in meters
    altitude: float = 828000.0

    # Antenna pattern of the sensor
    # Possible values: "ITU-R RS.1813", "ITU-R RS.1861 9a", "ITU-R RS.1861 9b",
    # "ITU-R RS.1861 9c", "ITU-R RS.2043", "OMNI"
    # TODO: check `x` and `y`:
    # @important: for EESS Passive, antenna pattern is from Recommendatio ITU-R RS.1813
    antenna_pattern: str = "ITU-R RS.1813"  # Antenna radiation pattern

    # Antenna efficiency for pattern described in ITU-R RS.1813 [0-1]
    # Efficiency factor for the antenna, range from 0 to 1
    antenna_efficiency: float = 0.6

    # Antenna diameter for ITU-R RS.1813 [m]
    antenna_diameter: float = 2.2  # Diameter of the antenna in meters

    # Receive antenna gain - applicable for 9a, 9b and OMNI [dBi]
    antenna_gain: float = 52.0  # Gain of the antenna in dBi

    # Channel model, possible values are "FSPL" (free-space path loss), "P619"
    channel_model: str = "FSPL"  # Channel model to be used

    # Parameters for the P.619 propagation model
    #    earth_station_alt_m - altitude of IMT system (in meters)
    #    earth_station_lat_deg - latitude of IMT system (in degrees)
    #    earth_station_long_diff_deg - difference between longitudes of IMT and satellite system
    #      (positive if space-station is to the East of earth-station)
    #    season - season of the year.
    param_p619 = ParametersP619()
    earth_station_alt_m: float = 0.0
    earth_station_lat_deg: float = 0.0
    earth_station_long_diff_deg: float = 0.0
    season: str = "SUMMER"

    ########### Creates a statistical distribution of nadir angle###############
    ############## following variables nadir_angle_distribution#################
    # if distribution_enable = ON, nadir_angle will vary statistically#########
    # if distribution_enable = OFF, nadir_angle follow nadir_angle variable ###
    # distribution_type = UNIFORM
    # UNIFORM = UNIFORM distribution in nadir_angle
    # 			- nadir_angle_distribution = initial nadir angle, final nadir angle
    distribution_enable: bool = False
    distribution_type: str = "UNIFORM"
    nadir_angle_distribution: tuple = (18.5, 49.3)

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
        # self.param_p619.load_from_paramters(self)
        # # print("altitude, earth_station_alt_m", self.altitude, self.earth_station_alt_m)
        # print([
        #     self.earth_station_alt_m, self.earth_station_lat_deg, self.earth_station_long_diff_deg
        # ])
        # Implement additional sanity checks for EESS specific parameters
        if not (0 <= self.antenna_efficiency or self.antenna_efficiency <= 1):
            raise ValueError("antenna_efficiency must be between 0 and 1")

        if self.antenna_pattern not in [
            "ITU-R RS.1813", "ITU-R RS.1861 9a",
            "ITU-R RS.1861 9b", "ITU-R RS.1861 9c",
            "ITU-R RS.2043", "OMNI",
        ]:
            raise ValueError(f"Invalid antenna_pattern: {self.antenna_pattern}")

        if self.antenna_pattern == "ITU-R RS.2043" and \
                (self.frequency <= 9000.0 or self.frequency >= 10999.0):
            raise ValueError(
                f"Frequency {self.frequency} MHz is not in the range for antenna pattern \"ITU-R RS.2043\"",
            )

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
