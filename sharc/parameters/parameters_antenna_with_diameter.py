from sharc.parameters.parameters_base import ParametersBase

from dataclasses import dataclass


@dataclass
class ParametersAntennaWithDiameter(ParametersBase):
    antenna_gain: float = None
    frequency: float = None
    diameter: float = None

    # ideally we would make this happen at the ParametersBase.load_subparameters, but
    # since parser is a bit hacky, this is created without acces to system freq and antenna gain
    # so validation needs to happen manually afterwards
    def validate(self, ctx):
        if None in [
            self.antenna_gain,
            self.diameter,
            self.frequency,
        ]:
            raise ValueError(f"{ctx} needs to have all its parameters set")

        if not isinstance(self.antenna_gain, int) and not isinstance(self.antenna_gain, float):
            raise ValueError(f"{ctx}.antenna_gain needs to be a number")

        if (not isinstance(self.diameter, int) and not isinstance(self.diameter, float)) or self.diameter <= 0:
            raise ValueError(f"{ctx}.diameter needs to be a positive number")

        if not isinstance(self.frequency, int) and not isinstance(self.frequency, float):
            raise ValueError(f"{ctx}.frequency needs to be a number")
