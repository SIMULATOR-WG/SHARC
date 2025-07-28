from sharc.parameters.parameters_antenna import ParametersAntenna

from sharc.antenna.antenna_fss_ss import AntennaFssSs
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.antenna.antenna_f699 import AntennaF699
from sharc.antenna.antenna_f1891 import AntennaF1891
from sharc.antenna.antenna_m1466 import AntennaM1466
from sharc.antenna.antenna_rs1813 import AntennaRS1813
from sharc.antenna.antenna_rs1861_9a import AntennaRS1861_9A
from sharc.antenna.antenna_rs1861_9b import AntennaRS1861_9B
from sharc.antenna.antenna_rs1861_9c import AntennaRS1861_9C
from sharc.antenna.antenna_rs2043 import AntennaRS2043
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_rra7_3 import AntennaReg_RR_A7_3
from sharc.antenna.antenna_modified_s465 import AntennaModifiedS465
from sharc.antenna.antenna_s580 import AntennaS580
from sharc.antenna.antenna_s672 import AntennaS672
from sharc.antenna.antenna_s1528 import AntennaS1528
from sharc.antenna.antenna_s1855 import AntennaS1855
from sharc.antenna.antenna_sa509 import AntennaSA509
from sharc.antenna.antenna_s1528 import AntennaS1528, AntennaS1528Leo, AntennaS1528Taylor
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt


class AntennaFactory():
    """Factory class for creating antenna instances based on pattern parameters."""

    @staticmethod
    def create_antenna(
        antenna_params: ParametersAntenna,
        azimuth: float,
        elevation: float,
    ):
        """Create and return an antenna instance based on the provided parameters, azimuth, and elevation."""
        match antenna_params.pattern:
            case "OMNI":
                return AntennaOmni(antenna_params.gain)
            case "ITU-R F.699":
                return AntennaF699(antenna_params.itu_r_f_699)
            case "ITU-R-S.1528-Taylor":
                return AntennaS1528Taylor(antenna_params.itu_r_s_1528)
            case "ITU-R-S.1528-LEO":
                return AntennaS1528Leo(antenna_params.itu_r_s_1528)
            case "ITU-R-S.1528-Section1.2":
                return AntennaS1528(antenna_params.itu_r_s_1528)
            case "ITU-R S.465":
                return AntennaS465(antenna_params.itu_r_s_465)
            case "ITU-R S.580":
                return AntennaS580(antenna_params.itu_r_s_580)
            case "MODIFIED ITU-R S.465":
                return AntennaS465(antenna_params.itu_r_s_465_modified)
            case "ITU-R S.1855":
                return AntennaS1855(antenna_params.itu_r_s_1855)
            case "ITU-R Reg. RR. Appendice 7 Annex 3":
                return AntennaReg_RR_A7_3(antenna_params.itu_reg_rr_a7_3)
            case "ARRAY":
                return AntennaBeamformingImt(
                    antenna_params.array.get_antenna_parameters(),
                    azimuth,
                    elevation
                )
            case _:
                raise ValueError(
                    f"Antenna factory does not support pattern {antenna_params.pattern}"
                )
