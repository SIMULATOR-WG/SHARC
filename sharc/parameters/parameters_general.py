# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sharc.sharc_definitions import SHARC_IMPLEMENTED_SYSTEMS
from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersGeneral(ParametersBase):
    """Dataclass containing the general parameters for the simulator
    """
    section_name: str = "general"
    num_snapshots: int = 10000
    imt_link: str = "DOWNLINK"
    system: str = "RAS"
    enable_cochannel: bool = False
    enable_adjacent_channel: bool = True
    seed: int = 101
    overwrite_output: bool = True
    output_dir: str = "output"
    output_dir_prefix: str = "output"

    def load_parameters_from_file(self, config_file: str):
        """Load the parameters from file an run a sanity check

        Parameters
        ----------
        file_name : str
            the path to the configuration file

        Raises
        ------
        ValueError
            if a parameter is not valid
        """
        super().load_parameters_from_file(config_file)

        # Now do the sanity check for some parameters
        if self.imt_link.upper() not in ["DOWNLINK", "UPLINK"]:
            raise ValueError(f"ParametersGeneral: \
                             Invalid value for parameter imt_link - {self.imt_link} \
                             Possible values are DOWNLINK and UPLINK")

        if self.system not in SHARC_IMPLEMENTED_SYSTEMS:
            raise ValueError(f"Invalid system name {self.system}")
