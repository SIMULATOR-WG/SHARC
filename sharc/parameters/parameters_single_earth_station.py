from dataclasses import dataclass, field
import typing

from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_antenna import ParametersAntenna
from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.parameters_p452 import ParametersP452
from sharc.parameters.parameters_hdfss import ParametersHDFSS


@dataclass
class ParametersSingleEarthStation(ParametersBase):
    """
    Defines parameters for passive Earth Exploration Satellite Service (EESS) sensors
    and their interaction with other services based on ITU recommendations.
    """
    section_name: str = "single_earth_station"
    nested_parameters_enabled: bool = True

    # for when other system needs it
    central_latitude: float = None
    central_altitude: float = None
    central_longitude: float = None

    # NOTE: when using P.619 it is suggested that polarization loss = 3dB is normal
    # also,
    # NOTE: Verification needed:
    # polarization mismatch between IMT BS and linear polarization = 3dB in earth P2P case?
    # = 0 is safer choice
    polarization_loss: float = 3.0

    # Sensor center frequency [MHz]
    frequency: float = None  # Center frequency of the sensor in MHz

    # Sensor bandwidth [MHz]
    bandwidth: float = None  # Bandwidth of the sensor in MHz

    # System receive noise temperature [K]
    noise_temperature: float = None

    # Adjacent channel selectivity [dB]
    adjacent_ch_selectivity: float = None

    # Adjacent channel selectivity [dB]
    adjacent_ch_emissions: float = None

    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density: float = None

    # Antenna pattern of the sensor
    antenna: ParametersAntenna = field(default_factory=ParametersAntenna)

    # Channel model, possible values are "FSPL" (free-space path loss), "P619"
    channel_model: typing.Literal[
        "FSPL", "P619",
        "P452",
    ] = "FSPL"  # Channel model to be used

    param_p619: ParametersP619 = field(default_factory=ParametersP619)
    # TODO: remove season from system parameter and put it as p619 parameter
    season: typing.Literal["WINTER", "SUMMER"] = "SUMMER"

    param_p452: ParametersP452 = field(default_factory=ParametersP452)

    param_hdfss: ParametersHDFSS = field(default_factory=ParametersHDFSS)

    @dataclass
    class EarthStationGeometry(ParametersBase):
        height: float = None

        @dataclass
        class FixedOrUniformDist(ParametersBase):
            __EXISTING_TYPES = ["UNIFORM_DIST", "FIXED"]
            type: typing.Literal["UNIFORM_DIST", "FIXED"] = None
            fixed: float = None

            @dataclass
            class UniformDistParams(ParametersBase):
                min: float = None
                max: float = None

                def validate(self, ctx):
                    if not isinstance(self.min, int) and not isinstance(self.min, float):
                        raise ValueError(
                            f"{ctx}.min parameter should be a number",
                        )
                    if not isinstance(self.max, int) and not isinstance(self.max, float):
                        raise ValueError(
                            f"{ctx}.max parameter should be a number",
                        )
                    if self.max < self.min:
                        raise ValueError(
                            f"{ctx}.max parameter should be greater than or equal to {ctx}.min",
                        )
            uniform_dist: UniformDistParams = field(
                default_factory=UniformDistParams,
            )

            def validate(self, ctx):
                if self.type not in self.__EXISTING_TYPES:
                    raise ValueError(
                        f"Invalid value for {ctx}.type. Should be one of {self.__EXISTING_TYPES}",
                    )

                match self.type:
                    case "UNIFORM_DIST":
                        self.uniform_dist.validate(f"{ctx}.uniform_dist")
                    case "FIXED":
                        if not isinstance(self.fixed, int) and not isinstance(self.fixed, float):
                            raise ValueError(f"{ctx}.fixed should be a number")
                    case _:
                        raise NotImplementedError(
                            f"Validation for {ctx}.type = {self.type} is not implemented",
                        )

        azimuth: FixedOrUniformDist = field(default_factory=FixedOrUniformDist)
        elevation: FixedOrUniformDist = field(
            default_factory=FixedOrUniformDist,
        )

        @dataclass
        class Location(ParametersBase):
            __EXISTING_TYPES = ["FIXED", "CELL", "NETWORK", "UNIFORM_DIST"]
            type: typing.Literal[
                "FIXED", "CELL",
                "NETWORK", "UNIFORM_DIST",
            ] = None

            @dataclass
            class LocationFixed(ParametersBase):
                x: float = None
                y: float = None

                def validate(self, ctx):
                    if not isinstance(self.x, int) and not isinstance(self.x, float):
                        raise ValueError(f"{ctx}.x needs to be a number")
                    if not isinstance(self.y, int) and not isinstance(self.y, float):
                        raise ValueError(f"{ctx}.y needs to be a number")

            @dataclass
            class LocationDistributed(ParametersBase):
                min_dist_to_bs: float = None

                def validate(self, ctx):
                    if not isinstance(self.min_dist_to_bs, int) and not isinstance(self.min_dist_to_bs, float):
                        raise ValueError(
                            f"{ctx}.min_dist_to_bs needs to be a number",
                        )

            @dataclass
            class LocationDistributedWithinCircle(ParametersBase):
                min_dist_to_center: float = None
                max_dist_to_center: float = None

                def validate(self, ctx):
                    if not isinstance(self.min_dist_to_center, int) and not isinstance(self.min_dist_to_center, float):
                        raise ValueError(
                            f"{ctx}.min_dist_to_center needs to be a number",
                        )
                    if not isinstance(self.max_dist_to_center, int) and not isinstance(self.max_dist_to_center, float):
                        raise ValueError(
                            f"{ctx}.max_dist_to_center needs to be a number",
                        )

            fixed: LocationFixed = field(default_factory=LocationFixed)
            cell: LocationDistributed = field(
                default_factory=LocationDistributed,
            )
            network: LocationDistributed = field(
                default_factory=LocationDistributed,
            )
            uniform_dist: LocationDistributedWithinCircle = field(
                default_factory=LocationDistributedWithinCircle,
            )

            def validate(self, ctx):
                match self.type:
                    case "FIXED":
                        self.fixed.validate(f"{ctx}.fixed")
                    case "CELL":
                        self.cell.validate(f"{ctx}.cell")
                    case "NETWORK":
                        self.network.validate(f"{ctx}.network")
                    case "UNIFORM_DIST":
                        self.uniform_dist.validate(f"{ctx}.uniform_dist")
                    case _:
                        raise NotImplementedError(
                            f"ParametersSingleEarthStation.Location.type = {self.type} has no validation implemented!",
                        )

        location: Location = field(default_factory=Location)

    geometry: EarthStationGeometry = field(
        default_factory=EarthStationGeometry,
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

        # this is needed because nested parameters
        # don't know/cannot access parents attributes
        self.antenna.set_external_parameters(
            frequency=self.frequency,
        )

        # this parameter is required in system get description
        self.antenna_pattern = self.antenna.pattern

        # this should be done by validating this parameters only if it is the selected system on the general section
        # TODO: make this better by changing the Parameters class itself
        should_validate = any(
            v is not None for v in [
                self.frequency, self.bandwidth,
            ]
        )

        if should_validate:
            self.validate(self.section_name)

    def validate(self, ctx="single_earth_station"):
        super().validate(ctx)

        if None in [self.frequency, self.bandwidth, self.channel_model]:
            raise ValueError(
                "ParametersSingleEarthStation required parameters are not all set",
            )

        if self.season not in ["WINTER", "SUMMER"]:
            raise ValueError(
                f"{ctx}.season needs to be either 'WINTER' or 'SUMMER'",
            )

        if self.channel_model not in ["FSPL", "P619", "P452", "TerrestrialSimple", "TVRO-URBAN", "TVRO-SUBURBAN"]:
            raise ValueError(
                f"{ctx}.channel_model" +
                "needs to be in ['FSPL', 'P619', 'P452', 'TerrestrialSimple', 'TVRO-URBAN', 'TVRO-SUBURBAN']",
            )
