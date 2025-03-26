from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_antenna_with_diameter import ParametersAntennaWithDiameter
from sharc.parameters.parameters_antenna_with_envelope_gain import ParametersAntennaWithEnvelopeGain

from dataclasses import dataclass, field
import typing


@dataclass
class ParametersAntenna(ParametersBase):
    # available antenna radiation patterns
    __SUPPORTED_ANTENNA_PATTERNS = [
        "OMNI", "ITU-R F.699", "ITU-R S.465", "ITU-R S.580", "MODIFIED ITU-R S.465", "ITU-R S.1855",
    ]

    # chosen antenna radiation pattern
    pattern: typing.Literal[
        "OMNI", "ITU-R F.699", "ITU-R S.465", "ITU-R S.580", "MODIFIED ITU-R S.465", "ITU-R S.1855",
    ] = None

    # antenna gain [dBi]
    gain: float = None

    itu_r_f_699: ParametersAntennaWithDiameter = field(
        default_factory=ParametersAntennaWithDiameter,
    )

    itu_r_s_465: ParametersAntennaWithDiameter = field(
        default_factory=ParametersAntennaWithDiameter,
    )

    itu_r_s_1855: ParametersAntennaWithDiameter = field(
        default_factory=ParametersAntennaWithDiameter,
    )

    itu_r_s_580: ParametersAntennaWithDiameter = field(
        default_factory=ParametersAntennaWithDiameter,
    )

    itu_r_s_465_modified: ParametersAntennaWithEnvelopeGain = field(
        default_factory=ParametersAntennaWithEnvelopeGain,
    )

    def set_external_parameters(self, *, frequency):
        attr_list = [
            a for a in dir(self) if not a.startswith('__') and isinstance(getattr(self, a), ParametersBase)
        ]

        for attr_name in attr_list:
            param = getattr(self, attr_name)

            if "frequency" in dir(param):
                param.frequency = frequency
            if "antenna_gain" in dir(param):
                param.antenna_gain = self.gain

    def load_parameters_from_file(self, config_file):
        # should not be loaded except as subparameter
        raise NotImplementedError()

    def validate(self, ctx):
        if None in [self.pattern, self.gain]:
            raise ValueError(
                f"Required parameters for Antenna were not met. At {ctx}",
            )

        if self.pattern not in self.__SUPPORTED_ANTENNA_PATTERNS:
            raise ValueError(
                f"Invalid ParametersAntenna.pattern. It should be one of: {self.__SUPPORTED_ANTENNA_PATTERNS}. At {ctx}",
            )

        match self.pattern:
            case "OMNI":
                pass
            case "ITU-R F.699":
                self.itu_r_f_699.validate(f"{ctx}.itu_r_f_699")
            case "ITU-R S.465":
                self.itu_r_s_465.validate(f"{ctx}.itu_r_s_465")
            case "ITU-R S.1855":
                self.itu_r_s_1855.validate(f"{ctx}.itu_r_s_1855")
            case "MODIFIED ITU-R S.465":
                self.itu_r_s_465_modified.validate(
                    f"{ctx}.itu_r_s_465_modified",
                )
            case "ITU-R S.580":
                self.itu_r_s_580.validate(f"{ctx}.itu_r_s_580")
            case _:
                raise NotImplementedError(
                    "ParametersAntenna.validate does not implement this antenna validation!",
                )
