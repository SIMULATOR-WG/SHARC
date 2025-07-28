from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersSingleBS(ParametersBase):
    intersite_distance: int = None
    cell_radius: int = None
    num_clusters: int = 1

    def load_subparameters(self, ctx: str, params: dict, quiet=True):
        super().load_subparameters(ctx, params, quiet)

        if self.intersite_distance is not None and self.cell_radius is not None:
            raise ValueError(f"{ctx}.intersite_distance and {ctx}.cell_radius should not be set at the same time")

        if self.intersite_distance is not None:
            self.cell_radius = self.intersite_distance * 2 / 3
        if self.cell_radius is not None:
            self.intersite_distance = self.cell_radius * 3 / 2

    def validate(self, ctx):
        if None in [self.intersite_distance, self.cell_radius]:
            raise ValueError(f"{ctx}.intersite_distance and cell_radius should be set.\
                                One of them through the parameters, the other inferred")

        if self.num_clusters not in [1, 2]:
            raise ValueError(f"{ctx}.num_clusters should either be 1 or 2")


if __name__ == "__main__":
    # Run validation for input parameters
    import os
    import pprint

    # Load default simulator parameters
    yaml_file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../../input/parameters.yaml")
    params = ParametersSingleBS()
    params.load_parameters_from_file(yaml_file_path, quiet=False)
    pprint.pprint(params)
