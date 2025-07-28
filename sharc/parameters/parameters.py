# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 19:35:52 2017

@author: edgar
"""

import sys
import os

from sharc.parameters.parameters_general import ParametersGeneral
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.parameters_eess_ss import ParametersEessSS
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_metsat_ss import ParametersMetSatSS
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_ras import ParametersRas
from sharc.parameters.parameters_single_earth_station import ParametersSingleEarthStation
from sharc.parameters.parameters_mss_ss import ParametersMssSs
from sharc.parameters.parameters_mss_d2d import ParametersMssD2d
from sharc.parameters.parameters_single_space_station import ParametersSingleSpaceStation


class Parameters(object):
    """
    Reads parameters from input file.
    """

    def __init__(self):
        self.file_name = None

        self.general = ParametersGeneral()
        self.imt = ParametersImt()
        self.eess_ss = ParametersEessSS()
        self.fs = ParametersFs()
        self.fss_ss = ParametersFssSs()
        self.fss_es = ParametersFssEs()
        self.haps = ParametersHaps()
        self.rns = ParametersRns()
        self.ras = ParametersRas()
        self.single_earth_station = ParametersSingleEarthStation()
        self.single_space_station = ParametersSingleSpaceStation()
        self.metsat_ss = ParametersMetSatSS()
        self.mss_ss = ParametersMssSs()
        self.mss_d2d = ParametersMssD2d()

    def set_file_name(self, file_name: str):
        """sets the configuration file name

        Parameters
        ----------
        file_name : str
            configuration file path
        """
        self.file_name = file_name

    def read_params(self):
        """Read the parameters from the config file
        """
        if not os.path.isfile(self.file_name):
            err_msg = f"PARAMETER ERROR [{self.__class__.__name__}]: \
                Could not find the configuration file {self.file_name}"
            sys.stderr.write(err_msg)
            sys.exit(1)

        #######################################################################
        # GENERAL
        #######################################################################
        self.general.load_parameters_from_file(self.file_name)

        #######################################################################
        # IMT
        #######################################################################
        self.imt.load_parameters_from_file(self.file_name)

        #######################################################################
        # FSS space station
        #######################################################################
        self.fss_ss.load_parameters_from_file(self.file_name)

        #######################################################################
        # FSS earth station
        #######################################################################
        self.fss_es.load_parameters_from_file(self.file_name)

        #######################################################################
        # Fixed wireless service
        #######################################################################
        self.fs.load_parameters_from_file(self.file_name)

        #######################################################################
        # HAPS (airbone) station
        #######################################################################
        self.haps.load_parameters_from_file(self.file_name)

        #######################################################################
        # RNS
        #######################################################################
        self.rns.load_parameters_from_file(self.file_name)

        #######################################################################
        # RAS station
        #######################################################################
        self.ras.load_parameters_from_file(self.file_name)

        #######################################################################
        # EESS passive
        #######################################################################
        self.eess_ss.load_parameters_from_file(self.file_name)

        #######################################################################
        # Single Earth Station
        #######################################################################
        self.single_earth_station.load_parameters_from_file(self.file_name)

        #######################################################################
        # MSS_SS
        #######################################################################
        self.mss_ss.load_parameters_from_file(self.file_name)

        #######################################################################
        # MSS_D2d
        #######################################################################
        self.mss_d2d.load_parameters_from_file(self.file_name)

        self.single_space_station.load_parameters_from_file(self.file_name)


if __name__ == "__main__":
    from pprint import pprint
    parameters = Parameters()
    param_sections = [
        a for a in dir(parameters) if not a.startswith('__') and not
        callable(getattr(parameters, a))
    ]
    print("\n#### Dumping default parameters:")
    for p in param_sections:
        print("\n")
        pprint(getattr(parameters, p))
