from dataclasses import dataclass, field
import typing

from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_antenna import ParametersAntenna
from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.parameters_p452 import ParametersP452
from sharc.parameters.parameters_hdfss import ParametersHDFSS


@dataclass
class ParametersSingleSpaceStation(ParametersBase):
    """
    Defines parameters for a single generic space station
    """
    section_name: str = "single_space_station"
    nested_parameters_enabled: bool = True
    is_space_to_earth: bool = True

    # Sensor center frequency [MHz]
    frequency: float = None  # Center frequency of the sensor in MHz

    # Sensor bandwidth [MHz]
    bandwidth: float = None  # Bandwidth of the sensor in MHz

    # System receive noise temperature [K]
    noise_temperature: float = None

    # Adjacent channel selectivity [dB]
    adjacent_ch_selectivity: float = None

    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density: float = None

    # Antenna pattern of the sensor
    antenna: ParametersAntenna = field(default_factory=ParametersAntenna)

    # Receiver polarization loss
    # e.g. could come from polarization mismatch or depolarization
    # check if IMT parameters don't come in values for single polarization
    # before adding loss here
    polarization_loss: float = 0.0

    # Channel model, possible values are "FSPL" (free-space path loss), "P619"
    channel_model: typing.Literal[
        "FSPL", "P619"
    ] = None  # Channel model to be used

    param_p619: ParametersP619 = field(default_factory=ParametersP619)
    # TODO: remove season from system parameter and put it as p619 parameter
    season: typing.Literal["WINTER", "SUMMER"] = "SUMMER"

    @dataclass
    class SpaceStationGeometry(ParametersBase):
        # NOTE: This does not directly translate to simulator 'height' param
        # default is GSO altitude
        altitude: float = 35786000.0
        es_altitude: float = 1200.0
        es_long_deg: float = -47.882778
        es_lat_deg: float = -15.793889

        @dataclass
        class PointingParam(ParametersBase):
            __EXISTING_TYPES = ["FIXED", "POINTING_AT_IMT"]
            type: typing.Literal["FIXED", "POINTING_AT_IMT"] = None
            fixed: float = None

            def validate(self, ctx):
                if self.type not in self.__EXISTING_TYPES:
                    raise ValueError(
                        f"Invalid value for {ctx}.type. Should be one of {self.__EXISTING_TYPES}",
                    )

                match self.type:
                    case "FIXED":
                        if not isinstance(self.fixed, int) and not isinstance(self.fixed, float):
                            raise ValueError(f"{ctx}.fixed should be a number")
                    case "POINTING_AT_IMT":
                        pass
                    case _:
                        raise NotImplementedError(
                            f"Validation for {ctx}.type = {self.type} is not implemented",
                        )

        azimuth: PointingParam = field(
            default_factory=lambda: ParametersSingleSpaceStation.SpaceStationGeometry.PointingParam(type="POINTING_AT_IMT"),
        )
        # default pointing directly downwards
        elevation: PointingParam = field(
            default_factory=lambda: ParametersSingleSpaceStation.SpaceStationGeometry.PointingParam(type="POINTING_AT_IMT"),
        )

        @dataclass
        class Location(ParametersBase):
            __EXISTING_TYPES = ["FIXED"]
            type: typing.Literal["FIXED"] = None

            @dataclass
            class LocationFixed(ParametersBase):
                # This should be the difference between SS lat/long and ES lat/long
                lat_deg: float = None
                long_deg: float = None

                def validate(self, ctx):
                    if not isinstance(self.lat_deg, int) and not isinstance(self.lat_deg, float):
                        raise ValueError(f"{ctx}.lat_deg needs to be a number")
                    if not isinstance(self.long_deg, int) and not isinstance(self.long_deg, float):
                        raise ValueError(f"{ctx}.long_deg needs to be a number")

            fixed: LocationFixed = field(default_factory=LocationFixed)

            def validate(self, ctx):
                match self.type:
                    case "FIXED":
                        self.fixed.validate(f"{ctx}.fixed")
                    case _:
                        raise NotImplementedError(
                            f"ParametersSingleSpaceStation.Location.type = {self.type} has no validation implemented!",
                        )

        location: Location = field(default_factory=Location)

    geometry: SpaceStationGeometry = field(
        default_factory=SpaceStationGeometry,
    )

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

        # this should be done by validating this parameters only if it is the selected system on the general section
        # TODO: make this better by changing the Parameters class itself
        should_validate = any(
            v is not None for v in [
                self.frequency, self.bandwidth,
            ]
        )

        if should_validate:
            self.validate(self.section_name)

    def propagate_parameters(self):
        if self.channel_model == "P619":
            if self.param_p619.earth_station_alt_m != ParametersP619.earth_station_alt_m:
                raise ValueError(
                    f"{self.section_name}.param_p619.earth_station_alt_m should not be set by hand."
                    "It is automatically set by other parameters in system"
                )
            if self.param_p619.earth_station_lat_deg != ParametersP619.earth_station_lat_deg:
                raise ValueError(
                    f"{self.section_name}.param_p619.earth_station_lat_deg should not be set by hand."
                    "It is automatically set by other parameters in system"
                )
            if self.param_p619.space_station_alt_m != ParametersP619.space_station_alt_m:
                raise ValueError(
                    f"{self.section_name}.param_p619.space_station_alt_m should not be set by hand."
                    "It is automatically set by other parameters in system"
                )
            if self.param_p619.earth_station_lat_deg != ParametersP619.earth_station_lat_deg:
                raise ValueError(
                    f"{self.section_name}.param_p619.earth_station_lat_deg should not be set by hand."
                    "It is automatically set by other parameters in system"
                )
        self.param_p619.space_station_alt_m = self.geometry.altitude
        self.param_p619.earth_station_alt_m = self.geometry.es_altitude
        self.param_p619.earth_station_lat_deg = self.geometry.es_lat_deg

        if self.geometry.location.type == "FIXED":
            self.param_p619.earth_station_long_diff_deg = self.geometry.location.fixed.long_deg - self.geometry.es_long_deg
        else:
            self.param_p619.earth_station_long_diff_deg = None

        # this is needed because nested parameters
        # don't know/cannot access parents attributes
        self.antenna.set_external_parameters(
            frequency=self.frequency,
        )

        # this parameter is required in system get description
        self.antenna_pattern = self.antenna.pattern

    def validate(self, ctx="single_space_station"):
        super().validate(ctx)

        if None in [self.frequency, self.bandwidth, self.channel_model, self.tx_power_density]:
            raise ValueError(
                "ParametersSingleSpaceStation required parameters are not all set",
            )

        if self.season not in ["WINTER", "SUMMER"]:
            raise ValueError(
                f"{ctx}.season needs to be either 'WINTER' or 'SUMMER'",
            )

        if self.channel_model not in ["FSPL", "P619"]:
            raise ValueError(
                f"{ctx}.channel_model" +
                "needs to be in ['FSPL', 'P619']",
            )
