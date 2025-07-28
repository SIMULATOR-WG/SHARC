# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 16:03:24 2017

@author: edgar
"""

import sys
import numpy.random as rnd
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.propagation.propagation_sat_simple import PropagationSatSimple
from sharc.propagation.propagation_ter_simple import PropagationTerSimple
from sharc.propagation.propagation_uma import PropagationUMa
from sharc.propagation.propagation_umi import PropagationUMi
from sharc.propagation.propagation_abg import PropagationABG
from sharc.propagation.propagation_clear_air_452 import PropagationClearAir
from sharc.propagation.propagation_tvro import PropagationTvro
from sharc.propagation.propagation_indoor import PropagationIndoor
from sharc.propagation.propagation_hdfss import PropagationHDFSS


class PropagationFactory(object):

    @staticmethod
    def create_propagation(
        channel_model: str,
        param: Parameters,
        param_system: ParametersBase,
        random_number_gen: rnd.RandomState,
    ) -> Propagation:
        """Creates a propagation model object

        Parameters
        ----------
        channel_model : str
            The channel model
        param : Parameters
            The simulation paramters.
        param_system : ParametersBase
            Specific system paramters. It can be either ParametersIMT or other system parameters.
        random_number_gen : rnd.RandomState
            Random number generator

        Returns
        -------
        Propagation
            Propagation object

        Raises
        ------
        ValueError
            Raises ValueError if the channel model is not implemented.
        """
        if channel_model == "FSPL":
            return PropagationFreeSpace(random_number_gen)
        elif channel_model == "ABG":
            return PropagationABG(random_number_gen)
        elif channel_model == "UMa":
            return PropagationUMa(random_number_gen)
        elif channel_model == "UMi":
            return PropagationUMi(random_number_gen, param.imt.los_adjustment_factor)
        elif channel_model == "SatelliteSimple":
            return PropagationSatSimple(random_number_gen)
        elif channel_model == "TerrestrialSimple":
            return PropagationTerSimple(random_number_gen)
        elif channel_model == "P619":
            if isinstance(param_system, ParametersImt):
                if param_system.topology.type != "NTN":
                    raise ValueError(
                        f"PropagationFactory: Channel model P.619 is invalid for topolgy {param.imt.topology.type}",
                    )
            else:
                # P.619 model is used only for space-to-earth links
                if param.imt.topology.type != "NTN" and not param_system.is_space_to_earth:
                    raise ValueError((
                        "PropagationFactory: Channel model P.619 is invalid"
                        f"for system {param.general.system} and IMT "
                        f"topology {param.imt.topology.type}"
                    ))
            return PropagationP619(
                random_number_gen=random_number_gen,
                space_station_alt_m=param_system.param_p619.space_station_alt_m,
                earth_station_alt_m=param_system.param_p619.earth_station_alt_m,
                earth_station_lat_deg=param_system.param_p619.earth_station_lat_deg,
                earth_station_long_diff_deg=param_system.param_p619.earth_station_lat_deg,
                season=param_system.season,
                mean_clutter_height=param_system.param_p619.mean_clutter_height,
                below_rooftop=param_system.param_p619.below_rooftop
            )
        elif channel_model == "P452":
            return PropagationClearAir(random_number_gen, param_system.param_p452)
        elif channel_model == "TVRO-URBAN":
            return PropagationTvro(random_number_gen, "URBAN")
        elif channel_model == "TVRO-SUBURBAN":
            return PropagationTvro(random_number_gen, "SUBURBAN")
        elif channel_model == "HDFSS":
            if param.general.system == "FSS_ES":
                # TODO: use param_hdfss in fss_es as well
                return PropagationHDFSS(param.fss_es, random_number_gen)
            else:
                return PropagationHDFSS(param_system.param_hdfss, random_number_gen)
        elif channel_model == "INDOOR":
            return PropagationIndoor(
                random_number_gen,
                param.imt.topology.indoor,
                param.imt.ue.k * param.imt.ue.k_m,
            )
        else:
            sys.stderr.write("ERROR\nInvalid channel_model: " + channel_model)
            sys.exit(1)
