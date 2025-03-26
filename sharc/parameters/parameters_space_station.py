from dataclasses import dataclass
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_p619 import ParametersP619
import math
from sharc.parameters.constants import EARTH_RADIUS


@dataclass
class ParametersSpaceStation(ParametersBase):
    """
    Defines parameters that should be used for orbiting Space Stations.
    TODO: use for FSS_SS in the future as well.
    """
    section_name: str = "space_station"
    is_space_to_earth: bool = True

    # Satellite center frequency [MHz]
    frequency: float = 0.0  # Center frequency of the satellite in MHz

    # Satellite bandwidth [MHz]
    bandwidth: float = 0.0  # Bandwidth of the satellite in MHz

    # Off-nadir pointing angle [deg]. Should only be set if elevation is not set
    nadir_angle: float = 0.0  # Angle in degrees away from the nadir point
    # Elevation angle [deg]. Should only be set if nadir_angle is not set
    elevation: float = 0.0

    # Satellite altitude [m]
    altitude: float = 0.0  # Altitude of the satellite above the Earth's surface in meters

    # satellite latitude [deg]
    lat_deg: float = 0.0

    # Antenna pattern of the satellite
    antenna_pattern: str | None = None  # Antenna radiation pattern

    # Antenna efficiency for pattern
    # Efficiency factor for the antenna, range from 0 to 1
    antenna_efficiency: float = 0.0

    # Antenna diameter
    antenna_diameter: float = 0.0  # Diameter of the antenna in meters

    # Receive antenna gain - applicable for 9a, 9b and OMNI [dBi]
    antenna_gain: float = 0.0  # Gain of the antenna in dBi

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

    # This parameter is also used by P619, but should not be set manually.
    # this should always == params.altitude
    space_station_alt_m: float = None

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

        if self.space_station_alt_m is not None:
            raise ValueError(
                "'space_station_alt_m' should not be set manually. It is always equal to 'altitude'",
            )

        if self.nadir_angle != 0 and self.elevation != 0:
            raise ValueError("'elevation' and 'nadir_angle' should not both be set at the same time. Choose either\
                             parameter to set")

        # Check channel model
        if self.channel_model not in ["FSPL", "P619"]:
            raise ValueError(
                "Invalid channel_model, must be either 'FSPL' or 'P619'",
            )

        if self.channel_model == "P619":
            # check necessary parameters for P619
            if None in [
                self.param_p619, self.earth_station_alt_m, self.earth_station_lat_deg, self.earth_station_long_diff_deg,
            ]:
                raise ValueError("When using P619 should set 'self.earth_station_alt_m', 'self.earth_station_lat_deg',\
                                 'self.earth_station_long_diff_deg' deafult or on parameter file and 'param_p619' \
                                 on child class default")
            # Check season
            if self.season not in ["SUMMER", "WINTER"]:
                raise ValueError(
                    "Invalid season, must be either 'SUMMER' or 'WINTER'",
                )

        self.set_derived_parameters()

    def set_derived_parameters(self):
        self.space_station_alt_m = self.altitude

        if self.param_p619:
            self.param_p619.load_from_paramters(self)

        if self.elevation != 0.0:
            # this relationship comes directly from law of sines
            self.nadir_angle = math.degrees(
                math.asin(
                    EARTH_RADIUS * math.sin(math.radians(self.elevation + 90)) /
                    (EARTH_RADIUS + self.altitude),
                ),
            )
        elif self.nadir_angle != 0.0:
            # this relationship comes directly from law of sines
            # can also be derived from incidence angle according to Rec. ITU-R RS.1861-0
            self.elevation = math.degrees(
                math.asin(
                    (EARTH_RADIUS + self.altitude) *
                    math.sin(math.radians(self.nadir_angle)) /
                    EARTH_RADIUS,
                ),
            ) - 90
