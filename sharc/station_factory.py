# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
"""

from warnings import warn
import numpy as np
import sys
import math

from sharc.support.enumerations import StationType
from sharc.parameters.parameters import Parameters
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_space_station import ParametersSpaceStation
from sharc.parameters.parameters_eess_ss import ParametersEessSS
from sharc.parameters.parameters_metsat_ss import ParametersMetSatSS
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_ras import ParametersRas
from sharc.parameters.parameters_single_earth_station import ParametersSingleEarthStation
from sharc.parameters.parameters_mss_ss import ParametersMssSs
from sharc.parameters.parameters_mss_d2d import ParametersMssD2d
from sharc.parameters.parameters_single_space_station import ParametersSingleSpaceStation
from sharc.parameters.constants import EARTH_RADIUS
from sharc.station_manager import StationManager
from sharc.mask.spectral_mask_imt import SpectralMaskImt
from sharc.antenna.antenna import Antenna
from sharc.antenna.antenna_factory import AntennaFactory
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
from sharc.antenna.antenna_multiple_transceiver import AntennaMultipleTransceiver
from sharc.topology.topology import Topology
from sharc.topology.topology_ntn import TopologyNTN
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.topology.topology_imt_mss_dc import TopologyImtMssDc
from sharc.mask.spectral_mask_3gpp import SpectralMask3Gpp
from sharc.mask.spectral_mask_mss import SpectralMaskMSS
from sharc.satellite.ngso.orbit_model import OrbitModel
from sharc.satellite.utils.sat_utils import calc_elevation, lla2ecef
from sharc.support.sharc_geom import rotate_angles_based_on_new_nadir, GeometryConverter

from sharc.parameters.constants import SPEED_OF_LIGHT


class StationFactory(object):

    @staticmethod
    def generate_imt_base_stations(
        param: ParametersImt,
        param_ant_bs: ParametersAntennaImt,
        topology: Topology,
        random_number_gen: np.random.RandomState,
    ):
        param_ant = param_ant_bs.get_antenna_parameters()
        num_bs = topology.num_base_stations
        imt_base_stations = StationManager(num_bs)
        imt_base_stations.station_type = StationType.IMT_BS
        if param.topology.type == "NTN":
            imt_base_stations.x = topology.space_station_x * np.ones(num_bs)
            imt_base_stations.y = topology.space_station_y * np.ones(num_bs)
            imt_base_stations.z = topology.space_station_z * np.ones(num_bs)
            imt_base_stations.height = imt_base_stations.z
            imt_base_stations.elevation = topology.elevation
            imt_base_stations.is_space_station = True
        elif param.topology.type == "MSS_DC":
            imt_base_stations.x = topology.space_station_x * np.ones(num_bs)
            imt_base_stations.y = topology.space_station_y * np.ones(num_bs)
            imt_base_stations.z = topology.space_station_z * np.ones(num_bs)
            imt_base_stations.height = topology.height
            imt_base_stations.elevation = topology.elevation
            imt_base_stations.is_space_station = True
        else:
            imt_base_stations.x = topology.x
            imt_base_stations.y = topology.y
            imt_base_stations.z = topology.z + param.bs.height
            imt_base_stations.elevation = -param_ant.downtilt * np.ones(num_bs)
            if param.topology.type == 'INDOOR':
                imt_base_stations.height = topology.height
            else:
                imt_base_stations.height = param.bs.height * np.ones(num_bs)

        imt_base_stations.azimuth = topology.azimuth
        imt_base_stations.active = random_number_gen.rand(
            num_bs,
        ) < param.bs.load_probability
        imt_base_stations.tx_power = param.bs.conducted_power * np.ones(num_bs)
        imt_base_stations.rx_power = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.rx_interference = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.ext_interference = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.total_interference = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )

        imt_base_stations.snr = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.sinr = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.sinr_ext = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )
        imt_base_stations.inr = dict(
            [(bs, -500 * np.ones(param.ue.k)) for bs in range(num_bs)],
        )

        imt_base_stations.antenna = np.empty(
            num_bs, dtype=Antenna,
        )

        for i in range(num_bs):
            imt_base_stations.antenna[i] = \
                AntennaFactory.create_antenna(
                    param.bs.antenna, imt_base_stations.azimuth[i],
                    imt_base_stations.elevation[i],)

        # imt_base_stations.antenna = [AntennaOmni(0) for bs in range(num_bs)]
        imt_base_stations.bandwidth = param.bandwidth * np.ones(num_bs)
        imt_base_stations.center_freq = param.frequency * np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs.noise_figure * \
            np.ones(num_bs)
        imt_base_stations.thermal_noise = -500 * np.ones(num_bs)

        if param.spectral_mask == "IMT-2020":
            imt_base_stations.spectral_mask = SpectralMaskImt(
                StationType.IMT_BS,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
                scenario=param.topology.type,
            )
        elif param.spectral_mask == "3GPP E-UTRA":
            imt_base_stations.spectral_mask = SpectralMask3Gpp(
                StationType.IMT_BS,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
            )

        if param.topology.type == 'MACROCELL':
            imt_base_stations.intersite_dist = param.topology.macrocell.intersite_distance
        elif param.topology.type == 'HOTSPOT':
            imt_base_stations.intersite_dist = param.topology.hotspot.intersite_distance

        return imt_base_stations

    @staticmethod
    def generate_imt_ue(
        param: ParametersImt,
        ue_param_ant: ParametersAntennaImt,
        topology: Topology,
        random_number_gen: np.random.RandomState,
    ) -> StationManager:

        if param.topology.type == "INDOOR":
            return StationFactory.generate_imt_ue_indoor(param, ue_param_ant, random_number_gen, topology)
        else:
            return StationFactory.generate_imt_ue_outdoor(param, ue_param_ant, random_number_gen, topology)

    @staticmethod
    def generate_ras_station(
        param: ParametersRas,
        random_number_gen: np.random.RandomState,
        topology: Topology,
    ) -> StationManager:
        return StationFactory.generate_single_earth_station(
            param, random_number_gen,
            StationType.RAS, topology
        )

    @staticmethod
    def generate_imt_ue_outdoor(
        param: ParametersImt,
        ue_param_ant: ParametersAntennaImt,
        random_number_gen: np.random.RandomState,
        topology: Topology,
    ) -> StationManager:
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue.k * param.ue.k_m

        num_ue = num_bs * num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE

        ue_x = list()
        ue_y = list()
        ue_z = list()

        imt_ue.height = param.ue.height * np.ones(num_ue)

        # TODO: Sanitaze the azimuth_range parameter
        azimuth_range = param.ue.azimuth_range
        if (not isinstance(azimuth_range, tuple)) or len(azimuth_range) != 2:
            raise ValueError("Invalid type or length for parameter azimuth_range")

        # Calculate UE pointing
        azimuth = (azimuth_range[1] - azimuth_range[0]) * \
            random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0]) * random_number_gen.random_sample(num_ue) + \
            elevation_range[0]

        if param.ue.distribution_type.upper() == "UNIFORM" or \
            param.ue.distribution_type.upper() == "CELL" or \
            param.ue.distribution_type.upper() == "UNIFORM_IN_CELL":

            central_cell = False
            deterministic_cell = False
            if param.ue.distribution_type.upper() != "UNIFORM":
                deterministic_cell = True
                if param.ue.distribution_type.upper() == "CELL":
                    central_cell = True
            else:
                if not (type(topology) is TopologyMacrocell):
                    sys.stderr.write(
                        "ERROR\nUniform UE distribution is currently supported only with Macrocell topology",
                    )
                    sys.exit(1)

            [ue_x, ue_y, ue_z, theta, distance] = StationFactory.get_random_position(
                num_ue, topology, random_number_gen,
                param.minimum_separation_distance_bs_ue,
                central_cell=central_cell,
                deterministic_cell=deterministic_cell,
            )
            psi = np.degrees(
                np.arctan((param.bs.height - param.ue.height) / distance),
            )

            imt_ue.azimuth = (azimuth + theta + np.pi / 2)
            imt_ue.elevation = elevation + psi

        elif param.ue.distribution_type.upper() == "ANGLE_AND_DISTANCE":
            # The Rayleigh and Normal distribution parameters (mean, scale and cutoff)
            # were agreed in TG 5/1 meeting (May 2017).

            if param.ue.distribution_distance.upper() == "RAYLEIGH":
                # For the distance between UE and BS, it is desired that 99% of UE's
                # are located inside the [soft] cell edge, i.e. Prob(d<d_edge) = 99%.
                # Since the distance is modeled by a random variable with Rayleigh
                # distribution, we use the quantile function to find that
                # sigma = distance/3.0345. So we always distibute UE's in order to meet
                # the requirement Prob(d<d_edge) = 99% for a given cell radius.
                radius_scale = topology.cell_radius / 3.0345
                radius = random_number_gen.rayleigh(radius_scale, num_ue)
                # clip the lower values to minimum_separation_distance_bs_ue
                radius = np.clip(radius, param.minimum_separation_distance_bs_ue, None)
            elif param.ue.distribution_distance.upper() == "UNIFORM":
                radius = (topology.cell_radius - param.minimum_separation_distance_bs_ue) * \
                    random_number_gen.random_sample(num_ue) - param.minimum_separation_distance_bs_ue
            else:
                sys.stderr.write(
                    "ERROR\nInvalid UE distance distribution: " + param.ue.distribution_distance,
                )
                sys.exit(1)

            if param.ue.distribution_azimuth.upper() == "NORMAL":
                # In case of the angles, we generate N times the number of UE's because
                # the angle cutoff will discard 5% of the terminals whose angle is
                # outside the angular sector defined by [-60, 60]. So, N = 1.4 seems to
                # be a safe choice.
                N = 1.4
                angle_scale = 30
                angle_mean = 0
                angle_n = random_number_gen.normal(
                    angle_mean, angle_scale, int(N * num_ue),
                )

                # FIXME: The angle cutoff ignores the lower range of azimuth_range
                angle_cutoff = np.max(azimuth_range)
                idx = np.where((angle_n < angle_cutoff) & (angle_n > -angle_cutoff))[0][:num_ue]
                angle = angle_n[idx]
            elif param.ue.distribution_azimuth.upper() == "UNIFORM":
                angle = (azimuth_range[1] - azimuth_range[0]) * random_number_gen.random_sample(num_ue) \
                    + azimuth_range[0]
            else:
                sys.stderr.write(
                    "ERROR\nInvalid UE azimuth distribution: " + param.ue.distribution_distance,
                )
                sys.exit(1)

            for bs in range(num_bs):
                idx = [
                    i for i in range(
                        bs * num_ue_per_bs, bs * num_ue_per_bs + num_ue_per_bs,
                    )
                ]
                # theta is the horizontal angle of the UE wrt the serving BS
                theta = topology.azimuth[bs] + angle[idx]
                # calculate UE position in x-y coordinates
                x = radius[idx] * np.cos(np.radians(theta))
                y = radius[idx] * np.sin(np.radians(theta))
                z = np.zeros_like(x)
                x, y, z = topology.transform_ue_xyz(
                    bs, x, y, z
                )
                ue_x.extend(x)
                ue_y.extend(y)
                ue_z.extend(z)

                # calculate UE azimuth wrt serving BS
                imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180) % 360

                # calculate elevation angle
                # psi is the vertical angle of the UE wrt the serving BS
                distance = np.sqrt(
                    (topology.x[bs] - x) ** 2 + (topology.y[bs] - y) ** 2,
                )
                psi = np.degrees(
                    np.arctan((param.bs.height - param.ue.height) / distance),
                )
                imt_ue.elevation[idx] = elevation[idx] + psi

        else:
            sys.stderr.write(
                "ERROR\nInvalid UE distribution type: " + param.ue.distribution_type,
            )
            sys.exit(1)

        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)
        imt_ue.z = np.array(ue_z) + param.ue.height

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.indoor = random_number_gen.random_sample(
            num_ue,
        ) <= (param.ue.indoor_percent / 100)
        imt_ue.rx_interference = -500 * np.ones(num_ue)
        imt_ue.ext_interference = -500 * np.ones(num_ue)

        # TODO: this piece of code works only for uplink
        par = ue_param_ant.get_antenna_parameters()
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaFactory.create_antenna(
                param.ue.antenna, imt_ue.azimuth[i],
                imt_ue.elevation[i],
            )

        # imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth * np.ones(num_ue)
        imt_ue.center_freq = param.frequency * np.ones(num_ue)
        imt_ue.noise_figure = param.ue.noise_figure * np.ones(num_ue)

        if param.spectral_mask == "IMT-2020":
            imt_ue.spectral_mask = SpectralMaskImt(
                StationType.IMT_UE,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
                scenario="OUTDOOR",
            )

        elif param.spectral_mask == "3GPP E-UTRA":
            imt_ue.spectral_mask = SpectralMask3Gpp(
                StationType.IMT_UE,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
            )

        imt_ue.spectral_mask.set_mask()

        if param.topology.type == 'MACROCELL':
            imt_ue.intersite_dist = param.topology.macrocell.intersite_distance
        elif param.topology.type == 'HOTSPOT':
            imt_ue.intersite_dist = param.topology.hotspot.intersite_distance

        return imt_ue

    @staticmethod
    def generate_imt_ue_indoor(
        param: ParametersImt,
        ue_param_ant: ParametersAntennaImt,
        random_number_gen: np.random.RandomState,
        topology: Topology,
    ) -> StationManager:
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue.k * param.ue.k_m
        num_ue = num_bs * num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE
        ue_x = list()
        ue_y = list()
        ue_z = list()

        # initially set all UE's as indoor
        imt_ue.indoor = np.ones(num_ue, dtype=bool)

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0]) * \
            random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        # azimuth = np.zeros(num_ue)
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0]) * \
            random_number_gen.random_sample(num_ue) + elevation_range[0]

        delta_x = (
            topology.b_w / math.sqrt(topology.ue_indoor_percent) - topology.b_w
        ) / 2
        delta_y = (
            topology.b_d / math.sqrt(topology.ue_indoor_percent) - topology.b_d
        ) / 2

        for bs in range(num_bs):
            idx = [
                i for i in range(
                    bs * num_ue_per_bs, bs * num_ue_per_bs + num_ue_per_bs,
                )
            ]
            # Right most cell of first floor
            if bs % topology.num_cells == 0 and bs < topology.total_bs_level:
                x_min = topology.x[bs] - topology.cell_radius - delta_x
                x_max = topology.x[bs] + topology.cell_radius
            # Left most cell of first floor
            elif bs % topology.num_cells == topology.num_cells - 1 and bs < topology.total_bs_level:
                x_min = topology.x[bs] - topology.cell_radius
                x_max = topology.x[bs] + topology.cell_radius + delta_x
            # Center cells and higher floors
            else:
                x_min = topology.x[bs] - topology.cell_radius
                x_max = topology.x[bs] + topology.cell_radius

            # First floor
            if bs < topology.total_bs_level:
                y_min = topology.y[bs] - topology.b_d / 2 - delta_y
                y_max = topology.y[bs] + topology.b_d / 2 + delta_y
            # Higher floors
            else:
                y_min = topology.y[bs] - topology.b_d / 2
                y_max = topology.y[bs] + topology.b_d / 2

            x = (x_max - x_min) * \
                random_number_gen.random_sample(num_ue_per_bs) + x_min
            y = (y_max - y_min) * \
                random_number_gen.random_sample(num_ue_per_bs) + y_min
            z = [
                topology.height[bs] - topology.b_h +
                param.ue.height for k in range(num_ue_per_bs)
            ]
            ue_x.extend(x)
            ue_y.extend(y)
            ue_z.extend(z)

            # theta is the horizontal angle of the UE wrt the serving BS
            theta = np.degrees(
                np.arctan2(
                    y - topology.y[bs], x - topology.x[bs],
                ),
            )
            # calculate UE azimuth wrt serving BS
            imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180) % 360

            # calculate elevation angle
            # psi is the vertical angle of the UE wrt the serving BS
            distance = np.sqrt(
                (topology.x[bs] - x)**2 + (topology.y[bs] - y)**2,
            )
            psi = np.degrees(
                np.arctan((param.bs.height - param.ue.height) / distance),
            )
            imt_ue.elevation[idx] = elevation[idx] + psi

            # check if UE is indoor
            if bs % topology.num_cells == 0:
                out = (x < topology.x[bs] - topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d / 2) | \
                      (y < topology.y[bs] - topology.b_d / 2)
            elif bs % topology.num_cells == topology.num_cells - 1:
                out = (x > topology.x[bs] + topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d / 2) | \
                      (y < topology.y[bs] - topology.b_d / 2)
            else:
                out = (y > topology.y[bs] + topology.b_d / 2) | \
                      (y < topology.y[bs] - topology.b_d / 2)
            imt_ue.indoor[idx] = ~ out

        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)
        imt_ue.height = np.array(ue_z)

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.rx_interference = -500 * np.ones(num_ue)
        imt_ue.ext_interference = -500 * np.ones(num_ue)

        # TODO: this piece of code works only for uplink
        par = ue_param_ant.get_antenna_parameters()
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(
                par, imt_ue.azimuth[i],
                imt_ue.elevation[i], ue_param_ant.subarray
            )

        # imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth * np.ones(num_ue)
        imt_ue.center_freq = param.frequency * np.ones(num_ue)
        imt_ue.noise_figure = param.ue.noise_figure * np.ones(num_ue)

        if param.spectral_mask == "IMT-2020":
            imt_ue.spectral_mask = SpectralMaskImt(
                StationType.IMT_UE,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
                scenario="INDOOR",
            )

        elif param.spectral_mask == "3GPP E-UTRA":
            imt_ue.spectral_mask = SpectralMask3Gpp(
                StationType.IMT_UE,
                param.frequency,
                param.bandwidth,
                param.spurious_emissions,
            )

        imt_ue.spectral_mask.set_mask()

        return imt_ue

    @staticmethod
    def generate_system(
        parameters: Parameters,
        topology: Topology,
        random_number_gen: np.random.RandomState,
        geometry_converter=GeometryConverter()
    ):
        if parameters.imt.topology.type == 'MACROCELL':
            intersite_dist = parameters.imt.topology.macrocell.intersite_distance
        elif parameters.imt.topology.type == 'HOTSPOT':
            intersite_dist = parameters.imt.topology.hotspot.intersite_distance

        if parameters.general.system == "METSAT_SS":
            return StationFactory.generate_metsat_ss(parameters.metsat_ss)
        elif parameters.general.system == "EESS_SS":
            return StationFactory.generate_eess_space_station(parameters.eess_ss)
        elif parameters.general.system == "FSS_ES":
            return StationFactory.generate_fss_earth_station(parameters.fss_es, random_number_gen, topology)
        elif parameters.general.system == "SINGLE_EARTH_STATION":
            return StationFactory.generate_single_earth_station(parameters.single_earth_station, random_number_gen,
                                                                StationType.SINGLE_EARTH_STATION, topology)
        elif parameters.general.system == "SINGLE_SPACE_STATION":
            return StationFactory.generate_single_space_station(parameters.single_space_station)
        elif parameters.general.system == "RAS":
            return StationFactory.generate_ras_station(parameters.ras, random_number_gen, topology)
        elif parameters.general.system == "FSS_SS":
            return StationFactory.generate_fss_space_station(parameters.fss_ss)
        elif parameters.general.system == "FS":
            return StationFactory.generate_fs_station(parameters.fs)
        elif parameters.general.system == "HAPS":
            return StationFactory.generate_haps(parameters.haps, intersite_dist, random_number_gen)
        elif parameters.general.system == "RNS":
            return StationFactory.generate_rns(parameters.rns, random_number_gen)
        elif parameters.general.system == "MSS_SS":
            return StationFactory.generate_mss_ss(parameters.mss_ss)
        elif parameters.general.system == "MSS_D2D":
            return StationFactory.generate_mss_d2d(parameters.mss_d2d, random_number_gen, geometry_converter)
        else:
            sys.stderr.write(
                "ERROR\nInvalid system: " +
                parameters.general.system,
            )
            sys.exit(1)

    @staticmethod
    def generate_single_space_station(param: ParametersSingleSpaceStation, simplify_dist_to_y=True):
        """
        Creates a single satellite based on parameters.
        In case simplify_dist_to_y == True (default) satellite will be only on y axis
        """
        space_station = StationManager(1)
        space_station.station_type = StationType.SINGLE_SPACE_STATION
        space_station.is_space_station = True

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km = (EARTH_RADIUS + param.geometry.altitude) / 1000
        dist_imt_centre_earth_km = (
            EARTH_RADIUS + param.geometry.es_altitude
        ) / 1000

        # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
        sat_lat_rad = param.geometry.location.fixed.lat_deg * np.pi / 180.
        imt_long_diff_rad = (param.geometry.location.fixed.long_deg - param.geometry.es_long_deg) * np.pi / 180.
        x1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
        z1 = dist_sat_centre_earth_km * np.sin(sat_lat_rad)

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = param.geometry.es_lat_deg * np.pi / 180.
        space_station.x = np.array(
            [x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)],
        ) * 1000
        space_station.y = np.array([y1]) * 1000
        space_station.height = np.array([
            (
                z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad) -
                dist_imt_centre_earth_km
            ) * 1000,
        ])

        # putting on y axis
        if simplify_dist_to_y:
            space_station.y = np.sqrt(space_station.x * space_station.x + space_station.y * space_station.y)
            space_station.x = np.zeros_like(space_station.x)

        if param.geometry.azimuth.type == "POINTING_AT_IMT":
            space_station.azimuth = np.rad2deg(np.arctan2(-space_station.y, -space_station.x))
        elif param.geometry.azimuth.type == "FIXED":
            space_station.azimuth = param.geometry.azimuth.fixed
        else:
            raise ValueError(f"Did not recognize azimuth type of {param.geometry.azimuth.type}")

        if param.geometry.azimuth.type == "POINTING_AT_IMT":
            gnd_elev = np.rad2deg(np.arctan2(
                space_station.height,
                np.sqrt(space_station.y * space_station.y + space_station.x * space_station.x)
            ))
            space_station.elevation = -gnd_elev
        elif param.geometry.azimuth.type == "FIXED":
            space_station.elevation = param.geometry.elevation.fixed
        else:
            raise ValueError(f"Did not recognize elevation type of {param.geometry.elevation.type}")

        space_station.active = np.array([True])
        space_station.tx_power = np.array(
            [param.tx_power_density + 10 *
                math.log10(param.bandwidth * 1e6) + 30],
        )
        space_station.rx_interference = -500

        space_station.antenna = np.array([
            AntennaFactory.create_antenna(param.antenna, space_station.azimuth[0],
                                          space_station.elevation[0])
        ])

        space_station.z = space_station.height
        space_station.bandwidth = param.bandwidth
        space_station.noise_temperature = param.noise_temperature
        space_station.thermal_noise = -500
        space_station.total_interference = -500

        return space_station

    @staticmethod
    def generate_fss_space_station(param: ParametersFssSs):
        fss_space_station = StationManager(1)
        fss_space_station.station_type = StationType.FSS_SS
        fss_space_station.is_space_station = True

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km = (EARTH_RADIUS + param.altitude) / 1000
        dist_imt_centre_earth_km = (
            EARTH_RADIUS + param.earth_station_alt_m
        ) / 1000

        # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
        sat_lat_rad = param.lat_deg * np.pi / 180.
        imt_long_diff_rad = param.earth_station_long_diff_deg * np.pi / 180.
        x1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
        z1 = dist_sat_centre_earth_km * np.sin(sat_lat_rad)

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = param.earth_station_lat_deg * np.pi / 180.
        fss_space_station.x = np.array(
            [x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)],
        ) * 1000
        fss_space_station.y = np.array([y1]) * 1000
        fss_space_station.height = np.array([
            (
                z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad) -
                dist_imt_centre_earth_km
            ) * 1000,
        ])
        fss_space_station.z = fss_space_station.height

        fss_space_station.azimuth = np.array([param.azimuth])
        fss_space_station.elevation = np.array([param.elevation])

        fss_space_station.active = np.array([True])
        fss_space_station.tx_power = np.array(
            [param.tx_power_density + 10 *
                math.log10(param.bandwidth * 1e6) + 30],
        )
        fss_space_station.rx_interference = np.array([-500])

        if param.antenna_pattern == "OMNI":
            fss_space_station.antenna = np.array(
                [AntennaOmni(param.antenna_gain)],
            )
        elif param.antenna_pattern == "ITU-R S.672":
            fss_space_station.antenna = np.array([AntennaS672(param)])
        elif param.antenna_pattern == "ITU-R S.1528":
            fss_space_station.antenna = np.array([AntennaS1528(param.antenna_s1528)])
        elif param.antenna_pattern == "FSS_SS":
            fss_space_station.antenna = np.array([AntennaFssSs(param)])
        else:
            sys.stderr.write(
                "ERROR\nInvalid FSS SS antenna pattern: " + param.antenna_pattern,
            )
            sys.exit(1)

        fss_space_station.bandwidth = np.array([param.bandwidth])
        fss_space_station.noise_temperature = np.array([param.noise_temperature])
        fss_space_station.thermal_noise = np.array([-500])
        fss_space_station.total_interference = np.array([-500])

        return fss_space_station

    @staticmethod
    def generate_fss_earth_station(param: ParametersFssEs, random_number_gen: np.random.RandomState, *args):
        """
        @deprecated

        Since this creates a Single Earth Station, you should use StationFactory.generate_single_earth_station instead.
        This will be deleted in the future.
        ----------------------------------
        Generates FSS Earth Station.

        Arguments:
            param: ParametersFssEs
            random_number_gen: np.random.RandomState
            topology (optional): Topology
        """
        warn(
            "This is deprecated, use StationFactory.generate_single_earth_station() instead; date=2024-10-11",
            DeprecationWarning, stacklevel=2,
        )

        if len(args):
            topology = args[0]

        fss_earth_station = StationManager(1)
        fss_earth_station.station_type = StationType.FSS_ES

        if param.location.upper() == "FIXED":
            fss_earth_station.x = np.array([param.x])
            fss_earth_station.y = np.array([param.y])
        elif param.location.upper() == "CELL":
            x, y, _, _ = StationFactory.get_random_position(
                1, topology, random_number_gen,
                param.min_dist_to_bs, True,
            )
            fss_earth_station.x = np.array(x)
            fss_earth_station.y = np.array(y)
        elif param.location.upper() == "NETWORK":
            x, y, _, _ = StationFactory.get_random_position(
                1, topology, random_number_gen,
                param.min_dist_to_bs, False,
            )
            fss_earth_station.x = np.array(x)
            fss_earth_station.y = np.array(y)
        elif param.location.upper() == "UNIFORM_DIST":
            # FSS ES is randomly (uniform) created inside a circle of radius
            # equal to param.max_dist_to_bs
            if param.min_dist_to_bs < 0:
                sys.stderr.write(
                    "ERROR\nInvalid minimum distance from FSS ES to BS: {}".format(
                        param.min_dist_to_bs,
                    ),
                )
                sys.exit(1)
            while (True):
                dist_x = random_number_gen.uniform(
                    -param.max_dist_to_bs, param.max_dist_to_bs,
                )
                dist_y = random_number_gen.uniform(
                    -param.max_dist_to_bs, param.max_dist_to_bs,
                )
                radius = np.sqrt(dist_x**2 + dist_y**2)
                if (radius > param.min_dist_to_bs) & (radius < param.max_dist_to_bs):
                    break
            fss_earth_station.x[0] = dist_x
            fss_earth_station.y[0] = dist_y
        else:
            sys.stderr.write(
                "ERROR\nFSS-ES location type {} not supported".format(
                    param.location),
            )
            sys.exit(1)

        fss_earth_station.z = np.array([param.height])
        fss_earth_station.height = np.array([param.height])

        if param.azimuth.upper() == "RANDOM":
            fss_earth_station.azimuth = np.array([random_number_gen.uniform(-180., 180.)])
        else:
            fss_earth_station.azimuth = np.array([float(param.azimuth)])

        elevation = random_number_gen.uniform(
            param.elevation_min, param.elevation_max,
        )
        fss_earth_station.elevation = np.array([elevation])

        fss_earth_station.active = np.array([True])
        fss_earth_station.tx_power = np.array(
            [param.tx_power_density + 10 *
                math.log10(param.bandwidth * 1e6) + 30],
        )
        fss_earth_station.rx_interference = np.array([-500])

        if param.antenna_pattern.upper() == "OMNI":
            fss_earth_station.antenna = np.array(
                [AntennaOmni(param.antenna_gain)],
            )
        elif param.antenna_pattern.upper() == "ITU-R S.1855":
            fss_earth_station.antenna = np.array([AntennaS1855(param)])
        elif param.antenna_pattern.upper() == "ITU-R S.465":
            fss_earth_station.antenna = np.array([AntennaS465(param)])
        elif param.antenna_pattern.upper() == "MODIFIED ITU-R S.465":
            fss_earth_station.antenna = np.array([AntennaModifiedS465(param)])
        elif param.antenna_pattern.upper() == "ITU-R S.580":
            fss_earth_station.antenna = np.array([AntennaS580(param)])
        else:
            sys.stderr.write(
                "ERROR\nInvalid FSS ES antenna pattern: " + param.antenna_pattern,
            )
            sys.exit(1)

        fss_earth_station.noise_temperature = np.array([param.noise_temperature])
        fss_earth_station.bandwidth = np.array([param.bandwidth])
        fss_earth_station.noise_temperature = np.array([param.noise_temperature])
        fss_earth_station.thermal_noise = np.array([-500])
        fss_earth_station.total_interference = np.array([-500])
        fss_earth_station.rx_interference = np.array([-500])
        fss_earth_station.rx_power = np.array([-500])

        return fss_earth_station

    @staticmethod
    def generate_single_earth_station(
        param: ParametersSingleEarthStation, random_number_gen: np.random.RandomState,
        station_type=StationType.SINGLE_EARTH_STATION, topology=None,
    ):
        """
        Generates a Single Earth Station.

        Arguments:
            param: ParametersSingleEarthStation
            random_number_gen: np.random.RandomState
            topology (optional): Topology
        """
        single_earth_station = StationManager(1)
        single_earth_station.station_type = StationType.SINGLE_EARTH_STATION

        match param.geometry.location.type:
            case "FIXED":
                single_earth_station.x = np.array(
                    [param.geometry.location.fixed.x],
                )
                single_earth_station.y = np.array(
                    [param.geometry.location.fixed.y],
                )
            case "CELL":
                x, y, _, _ = StationFactory.get_random_position(
                    1, topology, random_number_gen,
                    param.geometry.location.cell.min_dist_to_bs, True,
                )
                single_earth_station.x = np.array(x)
                single_earth_station.y = np.array(y)
            case "NETWORK":
                x, y, _, _ = StationFactory.get_random_position(
                    1, topology, random_number_gen,
                    param.geometry.location.network.min_dist_to_bs, False,
                )
                single_earth_station.x = np.array(x)
                single_earth_station.y = np.array(y)
            case "UNIFORM_DIST":
                # ES is randomly (uniform) created inside a circle of radius
                # equal to param.max_dist_to_bs
                if param.geometry.location.uniform_dist.min_dist_to_bs < 0:
                    sys.stderr.write(
                        "ERROR\nInvalid minimum distance from Single ES to BS: {}".format(
                            param.geometry.location.uniform_dist.min_dist_to_bs,
                        ),
                    )
                    sys.exit(1)
                while (True):
                    dist_x = random_number_gen.uniform(
                        -param.geometry.location.uniform_dist.max_dist_to_bs,
                        param.geometry.location.uniform_dist.max_dist_to_bs,
                    )
                    dist_y = random_number_gen.uniform(
                        -param.geometry.location.uniform_dist.max_dist_to_bs,
                        param.geometry.location.uniform_dist.max_dist_to_bs,
                    )
                    radius = np.sqrt(dist_x**2 + dist_y**2)
                    if (radius > param.geometry.location.uniform_dist.min_dist_to_bs) & \
                       (radius < param.geometry.location.uniform_dist.max_dist_to_bs):
                        break
                single_earth_station.x[0] = dist_x
                single_earth_station.y[0] = dist_y
            case _:
                sys.stderr.write(
                    "ERROR\nSingle-ES location type {} not supported".format(
                        param.geometry.location.type),
                )
                sys.exit(1)

        single_earth_station.z = np.array([param.geometry.height])
        single_earth_station.height = np.array([param.geometry.height])

        if param.geometry.azimuth.type == "UNIFORM_DIST":
            if param.geometry.azimuth.uniform_dist.min < -180:
                sys.stderr.write(
                    "ERROR\nInvalid minimum azimuth: {} < -180".format(
                        param.geometry.azimuth.uniform_dist.min),
                )
                sys.exit(1)
            if param.geometry.azimuth.uniform_dist.max > 180:
                sys.stderr.write(
                    "ERROR\nInvalid maximum azimuth: {} > 180".format(
                        param.geometry.azimuth.uniform_dist.max,
                    ),
                )
                sys.exit(1)
            single_earth_station.azimuth = np.array([
                random_number_gen.uniform(
                    param.geometry.azimuth.uniform_dist.min, param.geometry.azimuth.uniform_dist.max,
                ),
            ])
        else:
            single_earth_station.azimuth = np.array(
                [param.geometry.azimuth.fixed],
            )

        if param.geometry.elevation.type == "UNIFORM_DIST":
            single_earth_station.elevation = np.array([
                random_number_gen.uniform(
                    param.geometry.elevation.uniform_dist.min, param.geometry.elevation.uniform_dist.max,
                ),
            ])
        else:
            single_earth_station.elevation = np.array(
                [param.geometry.elevation.fixed],
            )

        match param.antenna.pattern:
            case "OMNI":
                single_earth_station.antenna = np.array(
                    [AntennaOmni(param.antenna.gain)],
                )
            case "ITU-R S.465":
                single_earth_station.antenna = np.array(
                    [AntennaS465(param.antenna.itu_r_s_465)],
                )
            case "ITU-R Reg. RR. Appendice 7 Annex 3":
                single_earth_station.antenna = np.array(
                    [AntennaReg_RR_A7_3(param.antenna.itu_reg_rr_a7_3)],
                )
            case "ITU-R S.1855":
                single_earth_station.antenna = np.array(
                    [AntennaS1855(param.antenna.itu_r_s_1855)],
                )
            case "MODIFIED ITU-R S.465":
                single_earth_station.antenna = np.array(
                    [AntennaModifiedS465(param.antenna.itu_r_s_465_modified)],
                )
            case "ITU-R S.580":
                single_earth_station.antenna = np.array(
                    [AntennaS580(param.antenna.itu_r_s_580)],
                )
            case _:
                sys.stderr.write(
                    "ERROR\nInvalid FSS ES antenna pattern: " + param.antenna_pattern,
                )
                sys.exit(1)

        single_earth_station.active = np.array([True])
        single_earth_station.bandwidth = np.array([param.bandwidth])

        single_earth_station.tx_power = np.array(
            [param.tx_power_density + 10 *
                math.log10(param.bandwidth * 1e6) + 30],
        )

        single_earth_station.noise_temperature = np.array([param.noise_temperature])

        # TODO: check why this would not be set on the StationManager() constructor itself?
        single_earth_station.rx_interference = np.array([-500])
        single_earth_station.thermal_noise = np.array([-500])
        single_earth_station.total_interference = np.array([-500])

        return single_earth_station

    @staticmethod
    def generate_fs_station(param: ParametersFs):
        """
        @deprecated
        Since this creates a Single Earth Station, you should use StationFactory.generate_single_earth_station instead.
        This will be deleted in the future
        """
        warn(
            "This is deprecated, use StationFactory.generate_single_earth_station() instead; date=2024-10-11",
            DeprecationWarning, stacklevel=2,
        )

        fs_station = StationManager(1)
        fs_station.station_type = StationType.FS

        fs_station.x = np.array([param.x])
        fs_station.y = np.array([param.y])
        fs_station.z = np.array([param.height])
        fs_station.height = np.array([param.height])

        fs_station.azimuth = np.array([param.azimuth])
        fs_station.elevation = np.array([param.elevation])

        fs_station.active = np.array([True])
        fs_station.tx_power = np.array(
            [param.tx_power_density + 10 *
                math.log10(param.bandwidth * 1e6) + 30],
        )
        fs_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            fs_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R F.699":
            fs_station.antenna = np.array([AntennaF699(param)])
        else:
            sys.stderr.write(
                "ERROR\nInvalid FS antenna pattern: " + param.antenna_pattern,
            )
            sys.exit(1)

        fs_station.noise_temperature = param.noise_temperature
        fs_station.bandwidth = np.array([param.bandwidth])

        return fs_station

    @staticmethod
    def generate_haps(param: ParametersHaps, intersite_distance: int, random_number_gen: np.random.RandomState):
        num_haps = 1
        haps = StationManager(num_haps)
        haps.station_type = StationType.HAPS
        haps.is_space_station = True

#        d = intersite_distance
#        h = (d/3)*math.sqrt(3)/2
#        haps.x = np.array([0, 7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
#        haps.y = np.array([0, 9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
        haps.x = np.array([0])
        haps.y = np.array([0])
        haps.z = param.altitude * np.ones(num_haps)

        haps.height = haps.z

        elev_max = 68.19  # corresponds to 50 km radius and 20 km altitude
        haps.azimuth = 360 * random_number_gen.random_sample(num_haps)
        haps.elevation = ((270 + elev_max) - (270 - elev_max)) * random_number_gen.random_sample(num_haps) + \
                         (270 - elev_max)

        haps.active = np.ones(num_haps, dtype=bool)

        haps.antenna = np.empty(num_haps, dtype=Antenna)

        if param.antenna_pattern == "OMNI":
            for i in range(num_haps):
                haps.antenna[i] = AntennaOmni(param.antenna_gain)
        elif param.antenna_pattern == "ITU-R F.1891":
            for i in range(num_haps):
                haps.antenna[i] = AntennaF1891(param)
        else:
            sys.stderr.write(
                "ERROR\nInvalid HAPS (airbone) antenna pattern: " +
                param.antenna_pattern,
            )
            sys.exit(1)

        haps.bandwidth = np.array([param.bandwidth])

        return haps

    @staticmethod
    def generate_rns(param: ParametersRns, random_number_gen: np.random.RandomState):
        num_rns = 1
        rns = StationManager(num_rns)
        rns.station_type = StationType.RNS
        rns.is_space_station = True

        rns.x = np.array([param.x])
        rns.y = np.array([param.y])
        rns.z = np.array([param.altitude])
        rns.height = np.array([param.altitude])

        # minimum and maximum values for azimuth and elevation
        azimuth = np.array([-30, 30])
        elevation = np.array([-30, 5])

        rns.azimuth = 90 + (azimuth[1] - azimuth[0]) * \
            random_number_gen.random_sample(num_rns) + azimuth[0]
        rns.elevation = (elevation[1] - elevation[0]) * \
            random_number_gen.random_sample(num_rns) + elevation[0]

        rns.active = np.ones(num_rns, dtype=bool)

        if param.antenna_pattern == "OMNI":
            rns.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R M.1466":
            rns.antenna = np.array(
                [AntennaM1466(param.antenna_gain, rns.azimuth, rns.elevation)],
            )
        else:
            sys.stderr.write(
                "ERROR\nInvalid RNS antenna pattern: " + param.antenna_pattern,
            )
            sys.exit(1)

        rns.bandwidth = np.array([param.bandwidth])
        rns.noise_temperature = param.noise_temperature
        rns.thermal_noise = -500
        rns.total_interference = -500
        rns.rx_interference = -500

        return rns

    @staticmethod
    def generate_eess_space_station(param: ParametersEessSS):
        if param.distribution_enable:
            if param.distribution_type == "UNIFORM":
                param.nadir_angle = np.random.uniform(
                    param.nadir_angle_distribution[0],
                    param.nadir_angle_distribution[1],
                )

        return StationFactory.generate_space_station(param, StationType.EESS_SS)

    @staticmethod
    def generate_metsat_ss(param: ParametersMetSatSS):
        return StationFactory.generate_space_station(param, StationType.METSAT_SS)

    @staticmethod
    def generate_space_station(param: ParametersSpaceStation, station_type: StationType):
        # this method uses off-nadir angle and altitude to infer the entire geometry
        # TODO: make this usable on more space station cases (initially only works for metsat and eess)
        space_station = StationManager(1)
        space_station.is_space_station = True
        space_station.station_type = station_type
        # assert param.station_type is not None

        incidence_angle = math.degrees(
            math.asin(
                math.sin(math.radians(param.nadir_angle)) *
                (1 + (param.altitude / EARTH_RADIUS)),
            ),
        )

        # distance to field of view centre according to Rec. ITU-R RS.1861-0
        distance = EARTH_RADIUS * \
            math.sin(math.radians(incidence_angle - param.nadir_angle)) / \
            math.sin(math.radians(param.nadir_angle))

        # Elevation at ground (centre of the footprint)
        theta_grd_elev = 90 - incidence_angle

        space_station.x = np.array([0])
        space_station.y = np.array(
            [distance * math.cos(math.radians(theta_grd_elev))],
        )
        space_station.height = np.array(
            [distance * math.sin(math.radians(theta_grd_elev))],
        )
        space_station.z = space_station.height

        # Elevation and azimuth at sensor wrt centre of the footprint
        # It is assumed the sensor is at y-axis, hence azimuth is 270 deg
        space_station.azimuth = np.array([270])
        space_station.elevation = np.array([-theta_grd_elev])

        space_station.active = np.array([True])
        space_station.rx_interference = np.array([-500])

        if param.antenna_pattern == "OMNI":
            space_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R RS.1813":
            space_station.antenna = np.array([AntennaRS1813(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9a":
            space_station.antenna = np.array([AntennaRS1861_9A(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9b":
            space_station.antenna = np.array([AntennaRS1861_9B(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9c":
            space_station.antenna = np.array([AntennaRS1861_9C()])
        elif param.antenna_pattern == "ITU-R RS.2043":
            space_station.antenna = np.array([AntennaRS2043()])
        elif param.antenna_pattern == "ITU-R S.672":
            space_station.antenna = np.array([AntennaS672(param)])
        else:
            sys.stderr.write(
                "ERROR\nInvalid EESS PASSIVE antenna pattern: " + param.antenna_pattern,
            )
            sys.exit(1)

        space_station.bandwidth = np.array([param.bandwidth])
        # Noise temperature is not an input parameter for yet used systems.
        # It is included here to calculate the useless I/N values
        space_station.noise_temperature = np.array([500.0])

        return space_station

    @staticmethod
    def generate_mss_ss(param_mss: ParametersMssSs):
        # We borrow the TopologyNTN geometry as it's the same for MSS_SS
        ntn_topology = TopologyNTN(param_mss.intersite_distance,
                                   param_mss.cell_radius,
                                   param_mss.altitude,
                                   param_mss.azimuth,
                                   param_mss.elevation,
                                   param_mss.num_sectors)
        ntn_topology.calculate_coordinates()

        num_bs = ntn_topology.num_base_stations
        mss_ss = StationManager(n=num_bs)
        mss_ss.station_type = StationType.MSS_SS
        mss_ss.x = ntn_topology.space_station_x * np.ones(num_bs) + param_mss.x
        mss_ss.y = ntn_topology.space_station_y * np.ones(num_bs) + param_mss.y
        mss_ss.z = ntn_topology.space_station_z * np.ones(num_bs)
        mss_ss.height = ntn_topology.space_station_z * np.ones(num_bs)
        mss_ss.elevation = ntn_topology.elevation
        mss_ss.is_space_station = True
        mss_ss.azimuth = ntn_topology.azimuth
        mss_ss.active = np.ones(num_bs, dtype=int)
        mss_ss.tx_power = np.ones(num_bs, dtype=int) * param_mss.tx_power_density + 10 * np.log10(param_mss.bandwidth * 10**6)
        mss_ss.antenna = np.empty(num_bs, dtype=AntennaS1528Leo)

        for i in range(num_bs):
            if param_mss.antenna_pattern == "ITU-R-S.1528-LEO":
                mss_ss.antenna[i] = AntennaS1528Leo(param_mss.antenna_s1528)
            elif param_mss.antenna_pattern == "ITU-R-S.1528-Section1.2":
                mss_ss.antenna[i] = AntennaS1528(param_mss.antenna_s1528)
            elif param_mss.antenna_pattern == "ITU-R-S.1528-Taylor":
                mss_ss.antenna[i] = AntennaS1528Taylor(param_mss.antenna_s1528)
            else:
                raise ValueError("generate_mss_ss: Invalid antenna type: {param_mss.antenna_pattern}")

        if param_mss.spectral_mask == "IMT-2020":
            mss_ss.spectral_mask = SpectralMaskImt(StationType.IMT_BS,
                                                   param_mss.frequency,
                                                   param_mss.bandwidth,
                                                   param_mss.spurious_emissions,
                                                   scenario="OUTDOOR")
        elif param_mss.spectral_mask == "3GPP E-UTRA":
            mss_ss.spectral_mask = SpectralMask3Gpp(StationType.IMT_BS,
                                                    param_mss.frequency,
                                                    param_mss.bandwidth,
                                                    param_mss.spurious_emissions,
                                                    scenario="OUTDOOR")
        elif params.spectral_mask == "MSS":
            mss_ss.spectral_mask = SpectralMaskMSS(params.frequency,
                                                   params.bandwidth,
                                                   params.spurious_emissions)
        else:
            raise ValueError(f"Invalid or not implemented spectral mask - {param_mss.spectral_mask}")
        mss_ss.spectral_mask.set_mask(param_mss.tx_power_density + 10 * np.log10(param_mss.bandwidth * 1e6))

        return mss_ss

    def generate_mss_d2d(
        params: ParametersMssD2d,
        random_number_gen: np.random.RandomState,
        geometry_converter: GeometryConverter,
        also_generate_inactive: bool = False
    ):
        """
        Generate the MSS D2D constellation with support for multiple orbits and base station visibility.

        Parameters
        ----------
        params : ParametersMssD2d
            Parameters for the MSS D2D system, including orbits and antenna configuration.
        random_number_gen : np.random.RandomState
            Random number generator for generating satellite positions.
        geometry_converter : GeometryConverter
            A converter that has already set a reference for coordinates transformation

        Returns
        -------
        StationManager
            A StationManager object containing satellite configurations and positions.
        """
        geometry_converter.validate()

        # Initialize the StationManager for the MSS D2D system
        mss_d2d_values = TopologyImtMssDc.get_coordinates(
            geometry_converter,
            params,
            random_number_gen,
            not also_generate_inactive
        )

        total_satellites = mss_d2d_values["num_satellites"]

        mss_d2d = StationManager(n=total_satellites)
        mss_d2d.station_type = StationType.MSS_D2D  # Set the station type to MSS D2D
        mss_d2d.is_space_station = True  # Indicate that the station is in space

        if params.spectral_mask == "IMT-2020":
            mss_d2d.spectral_mask = SpectralMaskImt(StationType.IMT_BS,
                                                    params.frequency,
                                                    params.bandwidth,
                                                    params.spurious_emissions,
                                                    scenario="OUTDOOR")
        elif params.spectral_mask == "3GPP E-UTRA":
            mss_d2d.spectral_mask = SpectralMask3Gpp(StationType.IMT_BS,
                                                     params.frequency,
                                                     params.bandwidth,
                                                     params.spurious_emissions,
                                                     scenario="OUTDOOR")
        elif params.spectral_mask == "MSS":
            mss_d2d.spectral_mask = SpectralMaskMSS(params.frequency,
                                                    params.bandwidth,
                                                    params.spurious_emissions)
        else:
            raise ValueError(f"Invalid or not implemented spectral mask - {params.spectral_mask}")
        mss_d2d.spectral_mask.set_mask(params.tx_power_density + 10 * np.log10(params.bandwidth * 1e6))

       # Configure satellite positions in the StationManager
        mss_d2d.x = mss_d2d_values["sat_x"]
        mss_d2d.y = mss_d2d_values["sat_y"]
        mss_d2d.z = mss_d2d_values["sat_z"]
        mss_d2d.elevation = mss_d2d_values["sat_antenna_elev"]
        mss_d2d.azimuth = mss_d2d_values["sat_antenna_azim"]
        mss_d2d.height = mss_d2d_values["sat_alt"]

        mss_d2d.active = np.zeros(total_satellites, dtype=bool)

        if mss_d2d_values["num_active_satellites"] != mss_d2d_values["num_satellites"]:
            mss_d2d.active[mss_d2d_values["active_satellites_idxs"]] = random_number_gen.uniform(
                size=len(mss_d2d_values["active_satellites_idxs"])
            ) < params.beams_load_factor
        else:
            # Set active satellite flags
            mss_d2d.active = random_number_gen.uniform(
                size=total_satellites
            ) < params.beams_load_factor

        # Initialize satellites antennas
        # we need to initialize them after coordinates transformation because of
        # repeated state (elevation and azimuth) inside multiple transceiver implementation
        mss_d2d.antenna = np.empty(total_satellites, dtype=AntennaS1528Leo)
        if params.antenna_pattern == "ITU-R-S.1528-LEO":
            antenna_pattern = AntennaS1528Leo(params.antenna_s1528)
        elif params.antenna_pattern == "ITU-R-S.1528-Section1.2":
            antenna_pattern = AntennaS1528(params.antenna_s1528)
        elif params.antenna_pattern == "ITU-R-S.1528-Taylor":
            antenna_pattern = AntennaS1528Taylor(params.antenna_s1528)
        else:
            raise ValueError("generate_mss_ss: Invalid antenna type: {param_mss.antenna_pattern}")

        for i in range(mss_d2d.num_stations):
            mss_d2d.antenna[i] = antenna_pattern

        return mss_d2d  # Return the configured StationManager

    @staticmethod
    def get_random_position(num_stas: int,
                            topology: Topology,
                            random_number_gen: np.random.RandomState,
                            min_dist_to_bs=0.,
                            central_cell=False,
                            deterministic_cell=False):
        """
        Generate UE random-possitions inside the topolgy area.

        Parameters
        ----------
        num_stas : int
            Number of UE stations
        topology : Topology
            The IMT topology object
        random_number_gen : np.random.RandomState
            Random number generator
        min_dist_to_bs : _type_, optional
            Minimum distance to the BS, by default 0.
        central_cell : bool, optional
            Whether the central cell in the cluster is used, by default False
        deterministic_cell : bool, optional
            Fix the cell to be used as anchor point, by default False

        Returns
        -------
        tuple
            x, y, z, azimuth and elevation angles.
        """
        hexagon_radius = topology.intersite_distance * 2 / 3

        x = np.array([])
        y = np.array([])
        z = np.array([])
        bs_x = -hexagon_radius
        bs_y = 0

        while len(x) < num_stas:
            num_stas_temp = num_stas - len(x)
            # generate UE uniformly in a triangle
            x_temp = random_number_gen.uniform(0, hexagon_radius * np.cos(np.pi / 6), num_stas_temp)
            y_temp = random_number_gen.uniform(0, hexagon_radius / 2, num_stas_temp)

            invert_index = np.arctan(y_temp / x_temp) > np.pi / 6
            y_temp[invert_index] = -(hexagon_radius / 2 - y_temp[invert_index])
            x_temp[invert_index] = (hexagon_radius * np.cos(np.pi / 6) - x_temp[invert_index])

            # randomly choose a hextant
            hextant = random_number_gen.random_integers(0, 5, num_stas_temp)
            hextant_angle = np.pi / 6 + np.pi / 3 * hextant

            old_x = x_temp
            x_temp = x_temp * np.cos(hextant_angle) - y_temp * np.sin(hextant_angle)
            y_temp = old_x * np.sin(hextant_angle) + y_temp * np.cos(hextant_angle)

            dist = np.sqrt((x_temp - bs_x) ** 2 + (y_temp - bs_y) ** 2)
            indices = dist > min_dist_to_bs

            x_temp = x_temp[indices]
            y_temp = y_temp[indices]

            x = np.append(x, x_temp)
            y = np.append(y, y_temp)

        x = x - bs_x
        y = y - bs_y

        # choose cells
        if central_cell:
            central_cell_indices = np.where((topology.x == 0) & (topology.y == 0))

            if not len(central_cell_indices[0]):
                sys.stderr.write("ERROR\nTopology does not have a central cell")
                sys.exit(1)

            cell = central_cell_indices[0][random_number_gen.random_integers(0, len(central_cell_indices[0]) - 1,
                                                                             num_stas)]
        elif deterministic_cell:
            num_bs = topology.num_base_stations
            stas_per_cell = num_stas / num_bs
            cell = np.repeat(np.arange(num_bs, dtype=int), stas_per_cell)

        else:  # random cells
            num_bs = topology.num_base_stations
            cell = random_number_gen.random_integers(0, num_bs - 1, num_stas)

        cell_x = topology.x[cell]
        cell_y = topology.y[cell]
        cell_z = topology.z[cell]

        # x = x + cell_x + hexagon_radius * np.cos(topology.azimuth[cell] * np.pi / 180)
        # y = y + cell_y + hexagon_radius * np.sin(topology.azimuth[cell] * np.pi / 180)
        old_x = x
        x = x * np.cos(np.radians(topology.azimuth[cell])) - y * np.sin(np.radians(topology.azimuth[cell]))
        y = old_x * np.sin(np.radians(topology.azimuth[cell])) + y * np.cos(np.radians(topology.azimuth[cell]))
        x = x + cell_x
        y = y + cell_y
        z = cell_z

        x = list(x)
        y = list(y)
        z = list(z)

        # calculate UE azimuth wrt serving BS
        if topology.is_space_station is False:
            theta = np.arctan2(y - cell_y, x - cell_x)

            # calculate elevation angle
            # psi is the vertical angle of the UE wrt the serving BS
            distance = np.sqrt((cell_x - x) ** 2 + (cell_y - y) ** 2)
        else:
            theta = np.arctan2(y - topology.space_station_y[cell], x - topology.space_station_x[cell])
            distance = np.sqrt((cell_x - x) ** 2 + (cell_y - y) ** 2 + (topology.bs.height)**2)

        return x, y, z, theta, distance


if __name__ == '__main__':
    rand_gen = np.random.RandomState(101)
    geometry_converter = GeometryConverter()

    # somente vou utilizar a translação que o satélite teoricamente sofreu:
    ref_lat = -14.1
    ref_long = -45.1
    ref_alt = 1200
    geometry_converter.set_reference(ref_lat, ref_long, ref_alt)
    from sharc.parameters.parameters_orbit import ParametersOrbit

    orbit = ParametersOrbit(
        n_planes=20,
        sats_per_plane=32,
        phasing_deg=3.9,
        long_asc_deg=18.0,
        inclination_deg=54.5,
        perigee_alt_km=525,
        apogee_alt_km=525
    )
    from sharc.parameters.imt.parameters_imt_mss_dc import ParametersImtMssDc
    params = ParametersImtMssDc(
        beam_radius=36516.0,
        num_beams=7,
        orbits=[orbit]
    )
    params.sat_is_active_if.conditions = ["MINIMUM_ELEVATION_FROM_ES"]
    params.sat_is_active_if.minimum_elevation_from_es = 5.0

    topology = TopologyImtMssDc(params, geometry_converter)

    topology.calculate_coordinates(rand_gen)

    topology.calculate_coordinates(rand_gen)

    parameters = ParametersImt()
    parameters.ue.k = 1
    parameters.ue.k_m = 1
    parameters.ue.azimuth_range = (-180, 180)
    parameters.ue.distribution_distance = "UNIFORM"
    parameters.ue.distribution_type = "ANGLE_AND_DISTANCE"
    parameters.ue.distribution_azimuth = "NORMAL"
    parameters.ue.height = 1.5
    parameters.ue.indoor_percent = 0
    parameters.bandwidth = 10
    parameters.frequency = 10
    parameters.ue.noise_figure = 0

    imt_ue = StationFactory.generate_imt_ue_outdoor(
        parameters,
        parameters.ue.antenna.array,
        rand_gen,
        topology
    )

    from sharc.satellite.scripts.plot_3d_param_file import plot_globe_with_borders
    fig = plot_globe_with_borders(True, geometry_converter)

    import plotly.graph_objects as go

    # fig.add_trace(go.Scatter3d(
    #     x=topology.x,
    #     y=topology.y,
    #     z=topology.z,
    #     mode='markers',
    #     marker=dict(size=3, color='green', opacity=0.8),
    #     showlegend=False
    # ))
    fig.add_trace(go.Scatter3d(
        x=topology.space_station_x,
        y=topology.space_station_y,
        z=topology.space_station_z,
        mode='markers',
        marker=dict(size=3, color='green', opacity=0.8),
        showlegend=False
    ))

    fig.add_trace(
        go.Scatter3d(
            x=imt_ue.x,
            y=imt_ue.y,
            z=imt_ue.z,
            mode='markers',
            marker=dict(size=1, color='red', opacity=1),
            showlegend=False
        )
    )

    # TODO: replace this with generate imt mss dc station
    st = StationManager(topology.num_base_stations)
    st.x = topology.space_station_x
    st.y = topology.space_station_y
    st.z = topology.space_station_z

    fig.add_trace(
        go.Scatter3d(
            x=[0],
            y=[0],
            z=[0],
            mode='markers',
            marker=dict(size=3, color='black', opacity=1),
            showlegend=False
        )
    )

    from sharc.support.sharc_geom import polar_to_cartesian
    # Plot beam boresight vectors
    boresight_length = 100 * 1e3  # Length of the boresight vectors for visualization
    boresight_x, boresight_y, boresight_z = polar_to_cartesian(
        boresight_length,
        imt_ue.azimuth,
        imt_ue.elevation
    )
    # Add arrow heads to the end of the boresight vectors
    for x, y, z, bx, by, bz in zip(imt_ue.x,
                                   imt_ue.y,
                                   imt_ue.z,
                                   boresight_x,
                                   boresight_y,
                                   boresight_z):
        fig.add_trace(go.Cone(
            x=[x + bx],
            y=[y + by],
            z=[z + bz],
            u=[bx],
            v=[by],
            w=[bz],
            colorscale=[[0, 'orange'], [1, 'orange']],
            sizemode='absolute',
            sizeref=2 * boresight_length / 5,
            showscale=False
        ))
    for x, y, z, bx, by, bz in zip(imt_ue.x,
                                   imt_ue.y,
                                   imt_ue.z,
                                   boresight_x,
                                   boresight_y,
                                   boresight_z):
        fig.add_trace(go.Scatter3d(
            x=[x, x + bx],
            y=[y, y + by],
            z=[z, z + bz],
            mode='lines',
            line=dict(color='orange', width=2),
            name='Boresight'
        ))
    # Suppress the legend for the boresight plot
    fig.update_traces(showlegend=False, selector=dict(name='Boresight'))

    # Maintain axis proportions
    fig.update_layout(scene_aspectmode='data')

    ref_x = imt_ue.x[11]
    ref_y = imt_ue.y[11]
    ref_z = imt_ue.z[11]
    range_scale = 1000

    range_scale = 5000
    ref_x = 0
    ref_y = 0
    ref_z = 0
    fig.update_layout(
        scene=dict(
            zaxis=dict(
                range=(-1e3 * range_scale + ref_z, 1e3 * range_scale + ref_z)
            ),
            yaxis=dict(
                range=(-1e3 * range_scale + ref_y, 1e3 * range_scale + ref_y)
            ),
            xaxis=dict(
                range=(-1e3 * range_scale + ref_x, 1e3 * range_scale + ref_x)
            ),
        )
    )

    # fig.tight_layout()
    fig.show()
