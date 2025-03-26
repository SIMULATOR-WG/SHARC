from sharc.parameters.parameters_base import ParametersBase

from dataclasses import dataclass


@dataclass
class ParametersAntennaWithEnvelopeGain(ParametersBase):
    antenna_gain: float = None
    envelope_gain: float = None

    # ideally we would make this happen at the ParametersBase.load_subparameters, but
    # since parser is a bit hacky, this is created without acces to system freq and antenna gain
    # so validation needs to happen manually afterwards
    def validate(self, ctx):
        if None in [
            self.antenna_gain,
            self.envelope_gain,
        ]:
            raise ValueError(f"{ctx} needs to have all its parameters set")

        if not isinstance(self.gain, int) and not isinstance(self.gain, float):
            raise ValueError(f"{ctx}.gain needs to be a number")

        if not isinstance(self.envelope_gain, int) and not isinstance(self.envelope_gain, float):
            raise ValueError(f"{ctx}.envelope_gain needs to be a number")
