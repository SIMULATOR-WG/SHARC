from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_antenna_with_diameter import ParametersAntennaWithDiameter
from sharc.parameters.parameters_antenna_with_envelope_gain import ParametersAntennaWithEnvelopeGain
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt

from dataclasses import dataclass, field
import typing


@dataclass
class ParametersAntenna(ParametersBase):
    # available antenna radiation patterns
    __SUPPORTED_ANTENNA_PATTERNS = [
        "OMNI", "ITU-R F.699", "ITU-R S.465", "ITU-R S.580", "MODIFIED ITU-R S.465", "ITU-R S.1855",
        "ITU-R Reg. RR. Appendice 7 Annex 3", "ARRAY", "ITU-R-S.1528-Taylor",
        "ITU-R-S.1528-Section1.2", "ITU-R-S.1528-LEO"
    ]

    # chosen antenna radiation pattern
    pattern: typing.Literal[
        "OMNI", "ITU-R F.699", "ITU-R S.465", "ITU-R S.580", "MODIFIED ITU-R S.465", "ITU-R S.1855",
        "ITU-R Reg. RR. Appendice 7 Annex 3", "ARRAY", "ITU-R-S.1528-Taylor",
        "ITU-R-S.1528-Section1.2", "ITU-R-S.1528-LEO"
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

    itu_reg_rr_a7_3: ParametersAntennaWithDiameter = field(
        default_factory=ParametersAntennaWithDiameter,
    )

    array: ParametersAntennaImt = field(default_factory=lambda: ParametersAntennaImt(downtilt=0.0))

    # TODO: maybe separate each different S.1528 parameter?
    itu_r_s_1528: ParametersAntennaS1528 = field(
        default_factory=ParametersAntennaS1528,
    )

    def set_external_parameters(self, **kwargs):
        attr_list = [
            a for a in dir(self) if not a.startswith('__') and isinstance(getattr(self, a), ParametersBase)
        ]

        for attr_name in attr_list:
            param = getattr(self, attr_name)

            for k, v in kwargs.items():
                if k in dir(param):
                    setattr(param, k, v)

            if "antenna_gain" in dir(param):
                param.antenna_gain = self.gain

    def load_parameters_from_file(self, config_file):
        # should not be loaded except as subparameter
        raise NotImplementedError()

    def validate(self, ctx):
        if None in [self.pattern]:
            raise ValueError(
                f"{ctx}.pattern should be set. Is None instead",
            )

        if self.pattern != "ARRAY" and self.gain is None:
            raise ValueError(
                f"{ctx}.gain should be set if not using array antenna.",
            )

        if self.pattern not in self.__SUPPORTED_ANTENNA_PATTERNS:
            raise ValueError(
                f"Invalid {ctx}.pattern. It should be one of: {self.__SUPPORTED_ANTENNA_PATTERNS}.",
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
            case "ITU-R Reg. RR. Appendice 7 Annex 3":
                if self.itu_reg_rr_a7_3.diameter is None:
                    # just hijacking validation since diameter is optional
                    self.itu_reg_rr_a7_3.diameter = 0
                self.itu_reg_rr_a7_3.validate(f"{ctx}.itu_reg_rr_a7_3")
            case "ARRAY":
                # TODO: validate here and make array non imt specific
                # self.array.validate(
                #     f"{ctx}.array",
                # )
                pass
            case "ITU-R-S.1528-Taylor":
                self.itu_r_s_1528.validate(f"{ctx}.itu_r_s_1528")
            case "ITU-R-S.1528-Section1.2":
                self.itu_r_s_1528.validate(f"{ctx}.itu_r_s_1528")
            case "ITU-R-S.1528-LEO":
                self.itu_r_s_1528.validate(f"{ctx}.itu_r_s_1528")
            case _:
                raise NotImplementedError(
                    "ParametersAntenna.validate does not implement this antenna validation!",
                )
