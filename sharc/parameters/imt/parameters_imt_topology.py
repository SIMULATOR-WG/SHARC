import typing
from dataclasses import field, dataclass

from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.imt.parameters_hotspot import ParametersHotspot
from sharc.parameters.imt.parameters_indoor import ParametersIndoor
from sharc.parameters.imt.parameters_macrocell import ParametersMacrocell
from sharc.parameters.imt.parameters_ntn import ParametersNTN
from sharc.parameters.imt.parameters_imt_mss_dc import ParametersImtMssDc
from sharc.parameters.imt.parameters_single_bs import ParametersSingleBS


@dataclass
class ParametersImtTopology(ParametersBase):
    nested_parameters_enabled = True

    type: typing.Literal[
        "MACROCELL", "HOTSPOT", "INDOOR", "SINGLE_BS", "NTN", "MSS_DC"
    ] = "MACROCELL"

    # these parameters are needed in case the other system requires coordinate
    # transformation
    central_latitude: float = None
    central_longitude: float = None
    central_altitude: float = None

    macrocell: ParametersMacrocell = field(default_factory=ParametersMacrocell)
    hotspot: ParametersHotspot = field(default_factory=ParametersHotspot)
    indoor: ParametersIndoor = field(default_factory=ParametersIndoor)
    single_bs: ParametersSingleBS = field(default_factory=ParametersSingleBS)
    ntn: ParametersNTN = field(default_factory=ParametersNTN)
    mss_dc: ParametersImtMssDc = field(default_factory=ParametersImtMssDc)

    def validate(self, ctx):
        match self.type:
            case "MACROCELL":
                self.macrocell.validate(f"{ctx}.macrocell")
            case "HOTSPOT":
                self.hotspot.validate(f"{ctx}.hotspot")
            case "INDOOR":
                self.indoor.validate(f"{ctx}.indoor")
            case "SINGLE_BS":
                self.single_bs.validate(f"{ctx}.single_bs")
            case "NTN":
                self.ntn.validate(f"{ctx}.ntn")
            case "MSS_DC":
                self.mss_dc.validate(f"{ctx}.mss_dc")
            case _:
                raise NotImplementedError(f"{ctx}.type == '{self.type}' may not be implemented")
