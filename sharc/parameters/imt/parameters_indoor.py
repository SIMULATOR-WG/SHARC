from dataclasses import dataclass
import typing

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersIndoor(ParametersBase):
    # Basic path loss model for indoor topology. Possible values:
    #       "FSPL" (free-space path loss),
    #       "INH_OFFICE" (3GPP Indoor Hotspot - Office)
    basic_path_loss: typing.Literal["INH_OFFICE", "FSPL"] = "INH_OFFICE"
    # Number of rows of buildings in the simulation scenario
    n_rows: int = 3
    # Number of colums of buildings in the simulation scenario
    n_colums: int = 2
    # Number of buildings containing IMT stations. Options:
    # 'ALL': all buildings contain IMT stations.
    # An integer representing the number of buildings.
    num_imt_buildings: typing.Union[typing.Literal["ALL"], int] = "ALL"
    # Street width (building separation) [m]
    street_width: int = 30
    # Intersite distance [m]
    intersite_distance: int = 40
    # Number of cells per floor
    num_cells: int = 3
    # Number of floors per building
    num_floors: int = 1
    # Percentage of indoor UE's [0, 1]
    ue_indoor_percent: int = .95
    # Building class: "TRADITIONAL" or "THERMALLY_EFFICIENT"
    building_class: typing.Literal[
        "TRADITIONAL", "THERMALLY_EFFICIENT"
    ] = "TRADITIONAL"
