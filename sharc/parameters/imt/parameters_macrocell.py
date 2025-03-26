from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersMacrocell(ParametersBase):
    intersite_distance: int = None
    wrap_around: bool = False
    num_clusters: int = 1

    def validate(self, ctx):
        if not isinstance(self.intersite_distance, int) and not isinstance(self.intersite_distance, float):
            raise ValueError(f"{ctx}.intersite_distance should be a number")

        if self.num_clusters not in [1, 7]:
            raise ValueError(f"{ctx}.num_clusters should be either 1 or 7")
