# NOTE: some parameters are without use in implementation as of this commit
# Some parameters should probably go to the BS geometry definition
# instead of being ntn only. Need to ensure the parameters make sense
# to other topologies as well, or validate that those parameters
# can't be set if another topology is chosen

from dataclasses import dataclass
import numpy as np

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersNTN(ParametersBase):
    """
    Simulation parameters for NTN network topology.
    """
    section_name: str = "ntn"

    # NTN Airborne Platform height (m)
    bs_height: float = None

    # NTN cell radius in network topology [m]
    cell_radius: float = 90000

    # NTN Intersite Distance (m). Intersite distance = Cell Radius * sqrt(3)
    # @important: for NTN, intersite distance means something different than normally,
    # since the BS's are placed at center of hexagons
    intersite_distance: float = None

    # BS azimuth
    # TODO: Put this elsewhere (in a bs.geometry for example) if needed by another model
    bs_azimuth: float = 45
    # BS elevation
    bs_elevation: float = 90

    # Number of sectors
    num_sectors: int = 7

    # TODO: implement the below parameters in the simulator. They are currently not used
    # Backoff Power [Layer 2] [dB]. Allowed: 7 sector topology - Layer 2
    bs_backoff_power: int = 3

    # NTN Antenna configuration
    bs_n_rows_layer1: int = 2
    bs_n_columns_layer1: int = 2
    bs_n_rows_layer2: int = 4
    bs_n_columns_layer2: int = 2

    def load_subparameters(self, ctx: str, params: dict, quiet=True):
        super().load_subparameters(ctx, params, quiet)

        if self.cell_radius is not None and self.intersite_distance is not None:
            raise ValueError(f"You cannot set both {ctx}.intersite_distance and {ctx}.cell_radius.")

        if self.cell_radius is not None:
            self.intersite_distance = self.cell_radius * np.sqrt(3)

        if self.intersite_distance is not None:
            self.cell_radius = self.intersite_distance / np.sqrt(3)

    def set_external_parameters(self, *, bs_height: float):
        """
            This method is used to "propagate" parameters from external context
            to the values required by ntn topology. It's not ideal, but it's done
            this way until we decide on a better way to model context.
        """
        self.bs_height = bs_height

    def validate(self, ctx: str):
        # Now do the sanity check for some parameters
        if self.num_sectors not in [1, 7, 19]:
            raise ValueError(
                f"ParametersNTN: Invalid number of sectors {self.num_sectors}",
            )

        if self.bs_height <= 0:
            raise ValueError(
                f"ParametersNTN: bs_height must be greater than 0, but is {self.bs_height}",
            )

        if self.cell_radius <= 0:
            raise ValueError(
                f"ParametersNTN: cell_radius must be greater than 0, but is {self.cell_radius}",
            )

        if self.intersite_distance <= 0:
            raise ValueError(
                f"ParametersNTN: intersite_distance must be greater than 0, but is {self.intersite_distance}",
            )

        if not isinstance(self.bs_backoff_power, int) or self.bs_backoff_power < 0:
            raise ValueError(
                f"ParametersNTN: bs_backoff_power must be a non-negative integer, but is {self.bs_backoff_power}",
            )

        if not np.all((0 <= self.bs_azimuth) & (self.bs_azimuth <= 360)):
            raise ValueError(
                "ParametersNTN: bs_azimuth values must be between 0 and 360 degrees",
            )

        if not np.all((0 <= self.bs_elevation) & (self.bs_elevation <= 90)):
            raise ValueError(
                "ParametersNTN: bs_elevation values must be between 0 and 90 degrees",
            )
