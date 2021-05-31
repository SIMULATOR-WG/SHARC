# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
@modified: Luciano Camilo Thu Marc 18 10:24:25 2021
"""

import numpy as np
import sys
import math

from sharc.support.enumerations import StationType
from sharc.parameters.parameters import Parameters
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_eess_passive import ParametersEessPassive
from sharc.parameters.parameters_ss_mleo import ParametersSsMLeo
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_arns import ParametersArns
from sharc.parameters.parameters_ras import ParametersRas
from sharc.parameters.parameters_hibs import ParametersHibs
from sharc.station_manager import StationManager
from sharc.mask.spectral_mask_imt import SpectralMaskImt
from sharc.antenna.antenna import Antenna
from sharc.antenna.antenna_fss_ss import AntennaFssSs
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.antenna.antenna_f699 import AntennaF699
from sharc.antenna.antenna_f1891 import AntennaF1891
from sharc.antenna.antenna_m1466 import AntennaM1466
from sharc.antenna.antenna_rs1813 import AntennaRS1813
from sharc.antenna.antenna_rs1861_9a import AntennaRS1861_9A
from sharc.antenna.antenna_rs1861_9b import AntennaRS1861_9B
from sharc.antenna.antenna_rs1861_9c import AntennaRS1861_9C
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_modified_s465 import AntennaModifiedS465
from sharc.antenna.antenna_s580 import AntennaS580
from sharc.antenna.antenna_s672 import AntennaS672
from sharc.antenna.antenna_s1528 import AntennaS1528
from sharc.antenna.antenna_s1528_leo import AntennaS1528_LEO
from sharc.antenna.antenna_s1855 import AntennaS1855
from sharc.antenna.antenna_sa509 import AntennaSA509
from sharc.antenna.antenna_omni_f1336 import AntennaOmniF1336
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.antenna.antenna_bessel import AntennaBessel
from sharc.antenna.antenna_f1245 import AntennaF1245
from sharc.antenna.antenna_cossecant_squared import AntennaCossecantSquared
from sharc.antenna.antenna_meteorological_radar_cosine_n1 import AntennaMeteorologicalRadarCosine1
from sharc.antenna.antenna_metereological_radar_uniform import AntennaMeteorologicalRadarUniform
from sharc.antenna.antenna_radar_phased_array import AntennaRadarPhasedArray
from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.topology.topology_hibs import TopologyHIBS
from sharc.mask.spectral_mask_3gpp import SpectralMask3Gpp


class StationFactory(object):

    @staticmethod
    def generate_imt_ue(param: ParametersImt,
                        param_ant: ParametersAntennaImt,
                        topology: Topology,
                        random_number_gen: np.random.RandomState) -> StationManager:

        if param.topology == "INDOOR":
            return StationFactory.generate_imt_ue_indoor(param, param_ant, random_number_gen, topology)
        else:
            return StationFactory.generate_imt_ue_outdoor(param, param_ant, random_number_gen, topology)

    @staticmethod
    def generate_imt_ue_outdoor(param: ParametersImt,
                                param_ant: ParametersAntennaImt,
                                random_number_gen: np.random.RandomState,
                                topology: Topology) -> StationManager:
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k * param.ue_k_m

        num_ue = num_bs * num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE

        ue_x = list()
        ue_y = list()

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0]) * random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0]) * random_number_gen.random_sample(num_ue) + \
                    elevation_range[0]

        if param.ue_distribution_type.upper() == "UNIFORM":

            if not (type(topology) is TopologyMacrocell or TopologyHIBS):
                sys.stderr.write("ERROR\nUniform UE distribution is currently supported only with Macrocell topology")
                sys.exit(1)

            [ue_x, ue_y, theta, distance] = StationFactory.get_random_position(num_ue, topology, random_number_gen,
                                                                               param.minimum_separation_distance_bs_ue,
                                                                               deterministic_cell=True)
            psi = np.degrees(np.arctan((param.bs_height - param.ue_height) / distance))

            imt_ue.azimuth = (azimuth + theta + np.pi / 2)
            imt_ue.elevation = elevation + psi

        elif param.ue_distribution_type.upper() == "ANGLE_AND_DISTANCE":
            # The Rayleigh and Normal distribution parameters (mean, scale and cutoff)
            # were agreed in TG 5/1 meeting (May 2017).

            if param.ue_distribution_distance.upper() == "RAYLEIGH":
                # For the distance between UE and BS, it is desired that 99% of UE's
                # are located inside the [soft] cell edge, i.e. Prob(d<d_edge) = 99%.
                # Since the distance is modeled by a random variable with Rayleigh
                # distribution, we use the quantile function to find that
                # sigma = distance/3.0345. So we always distibute UE's in order to meet
                # the requirement Prob(d<d_edge) = 99% for a given cell radius.
                radius_scale = topology.cell_radius / 3.0345
                radius = random_number_gen.rayleigh(radius_scale, num_ue)
            elif param.ue_distribution_distance.upper() == "UNIFORM":
                radius = topology.cell_radius * random_number_gen.random_sample(num_ue)
            else:
                sys.stderr.write("ERROR\nInvalid UE distance distribution: " + param.ue_distribution_distance)
                sys.exit(1)

            if param.ue_distribution_azimuth.upper() == "NORMAL":
                # In case of the angles, we generate N times the number of UE's because
                # the angle cutoff will discard 5% of the terminals whose angle is
                # outside the angular sector defined by [-60, 60]. So, N = 1.4 seems to
                # be a safe choice.
                N = 1.4
                angle_scale = 30
                angle_mean = 0
                angle_n = random_number_gen.normal(angle_mean, angle_scale, int(N * num_ue))

                angle_cutoff = 60
                idx = np.where((angle_n < angle_cutoff) & (angle_n > -angle_cutoff))[0][:num_ue]
                angle = angle_n[idx]
            elif param.ue_distribution_azimuth.upper() == "UNIFORM":
                azimuth_range = (-60, 60)
                angle = (azimuth_range[1] - azimuth_range[0]) * random_number_gen.random_sample(num_ue) \
                        + azimuth_range[0]
            else:
                sys.stderr.write("ERROR\nInvalid UE azimuth distribution: " + param.ue_distribution_distance)
                sys.exit(1)

            for bs in range(num_bs):
                idx = [i for i in range(bs * num_ue_per_bs, bs * num_ue_per_bs + num_ue_per_bs)]
                # theta is the horizontal angle of the UE wrt the serving BS
                theta = topology.azimuth[bs] + angle[idx]
                # calculate UE position in x-y coordinates
                x = topology.x[bs] + radius[idx] * np.cos(np.radians(theta))
                y = topology.y[bs] + radius[idx] * np.sin(np.radians(theta))
                ue_x.extend(x)
                ue_y.extend(y)

                # calculate UE azimuth wrt serving BS
                imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180) % 360

                # calculate elevation angle
                # psi is the vertical angle of the UE wrt the serving BS
                distance = np.sqrt((topology.x[bs] - x) ** 2 + (topology.y[bs] - y) ** 2)
                psi = np.degrees(np.arctan((param.bs_height - param.ue_height) / distance))
                imt_ue.elevation[idx] = elevation[idx] + psi
        else:
            sys.stderr.write("ERROR\nInvalid UE distribution type: " + param.ue_distribution_type)
            sys.exit(1)

        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.height = param.ue_height * np.ones(num_ue)
        imt_ue.indoor = random_number_gen.random_sample(num_ue) <= (param.ue_indoor_percent / 100)
        imt_ue.rx_interference = -500 * np.ones(num_ue)
        imt_ue.ext_interference = -500 * np.ones(num_ue)

        # TODO: this piece of code works only for uplink
        par = param_ant.get_antenna_parameters(StationType.IMT_UE)
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(par, imt_ue.azimuth[i],
                                                      imt_ue.elevation[i])

        # imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth * np.ones(num_ue)
        imt_ue.center_freq = param.frequency * np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure * np.ones(num_ue)

        if param.spectral_mask == "IMT-2020":
            imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,
                                                   param.frequency,
                                                   param.bandwidth,
                                                   param.spurious_emissions,
                                                   scenario="OUTDOOR")

        elif param.spectral_mask == "3GPP E-UTRA":
            imt_ue.spectral_mask = SpectralMask3Gpp(StationType.IMT_UE,
                                                    param.frequency,
                                                    param.bandwidth,
                                                    param.spurious_emissions)

        elif param.spectral_mask == "3GPP 36.104":
            imt_ue.spectral_mask = SpectralMask3Gpp(StationType.IMT_UE, param.frequency, param.bandwidth,
                                                    param.spurious_emissions)

        imt_ue.spectral_mask.set_mask()

        if param.topology == 'MACROCELL' or param.topology == 'HOTSPOT':
            imt_ue.intersite_dist = param.intersite_distance

        return imt_ue

    @staticmethod
    def generate_imt_ue_indoor(param: ParametersImt,
                               param_ant: ParametersAntennaImt,
                               random_number_gen: np.random.RandomState,
                               topology: Topology) -> StationManager:
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k * param.ue_k_m
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
        azimuth = (azimuth_range[1] - azimuth_range[0]) * random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        # azimuth = np.zeros(num_ue)
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0]) * random_number_gen.random_sample(num_ue) + \
                    elevation_range[0]

        delta_x = (topology.b_w / math.sqrt(topology.ue_indoor_percent) - topology.b_w) / 2
        delta_y = (topology.b_d / math.sqrt(topology.ue_indoor_percent) - topology.b_d) / 2

        for bs in range(num_bs):
            idx = [i for i in range(bs * num_ue_per_bs, bs * num_ue_per_bs + num_ue_per_bs)]
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

            x = (x_max - x_min) * random_number_gen.random_sample(num_ue_per_bs) + x_min
            y = (y_max - y_min) * random_number_gen.random_sample(num_ue_per_bs) + y_min
            z = [topology.height[bs] - topology.b_h + param.ue_height for k in range(num_ue_per_bs)]
            ue_x.extend(x)
            ue_y.extend(y)
            ue_z.extend(z)

            # theta is the horizontal angle of the UE wrt the serving BS
            theta = np.degrees(np.arctan2(y - topology.y[bs], x - topology.x[bs]))
            # calculate UE azimuth wrt serving BS
            imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180) % 360

            # calculate elevation angle
            # psi is the vertical angle of the UE wrt the serving BS
            distance = np.sqrt((topology.x[bs] - x) ** 2 + (topology.y[bs] - y) ** 2)
            psi = np.degrees(np.arctan((param.bs_height - param.ue_height) / distance))
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
        par = param_ant.get_antenna_parameters(StationType.IMT_UE)
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(par, imt_ue.azimuth[i],
                                                      imt_ue.elevation[i])

        # imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth * np.ones(num_ue)
        imt_ue.center_freq = param.frequency * np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure * np.ones(num_ue)

        if param.spectral_mask == "IMT-2020":
            imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,
                                                   param.frequency,
                                                   param.bandwidth,
                                                   param.spurious_emissions,
                                                   scenario="INDOOR")

        elif param.spectral_mask == "3GPP E-UTRA":
            imt_ue.spectral_mask = SpectralMask3Gpp(StationType.IMT_UE,
                                                    param.frequency,
                                                    param.bandwidth,
                                                    param.spurious_emissions)

        imt_ue.spectral_mask.set_mask()

        return imt_ue

    @staticmethod
    def generate_system(parameters: Parameters, topology: Topology, random_number_gen: np.random.RandomState):
        if parameters.general.system == "EESS_PASSIVE":
            return StationFactory.generate_eess_passive_sensor(parameters.eess_passive)
        if parameters.general.system == "FSS_ES":
            return StationFactory.generate_fss_earth_station(parameters.fss_es, random_number_gen, topology)
        elif parameters.general.system == "SS_MLEO":
            return StationFactory.generate_ss_mleo_space_station(parameters.ss_mleo)
        elif parameters.general.system == "FSS_SS":
            return StationFactory.generate_fss_space_station(parameters.fss_ss)
        elif parameters.general.system == "FS":
            return StationFactory.generate_fs_station(parameters.fs)
        elif parameters.general.system == "HAPS":
            return StationFactory.generate_haps(parameters.haps, parameters.imt.intersite_distance, random_number_gen)
        elif parameters.general.system == "RNS":
            return StationFactory.generate_rns(parameters.rns, random_number_gen)
        elif parameters.general.system == "ARNS":
            return StationFactory.generate_arns(parameters.arns, random_number_gen)
        elif parameters.general.system == "RAS":
            return StationFactory.generate_ras_station(parameters.ras)
        else:
            sys.stderr.write("ERROR\nInvalid system: " + parameters.general.system)
            sys.exit(1)

    @staticmethod
    def generate_ss_mleo_space_station(param: ParametersSsMLeo):
        ss_mleo_space_station = StationManager(1)
        ss_mleo_space_station.station_type = StationType.SS_MLEO

        # ss_mleo_space_station.x = np.array([param.x])
        # ss_mleo_space_station.y = np.array([param.y])
        # ss_mleo_space_station.height = np.array([param.height])

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km = (param.EARTH_RADIUS + param.altitude) / 1000
        dist_imt_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude) / 1000

        # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
        sat_lat_rad = param.lat_deg * np.pi / 180.
        imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
        x1 = dist_sat_centre_earth_km * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth_km * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
        z1 = dist_sat_centre_earth_km * np.sin(sat_lat_rad)

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = param.imt_lat_deg * np.pi / 180.
        ss_mleo_space_station.x = np.array([x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)]) * 1000
        ss_mleo_space_station.y = np.array([y1]) * 1000
        ss_mleo_space_station.height = np.array([(z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad)
                                                  - dist_imt_centre_earth_km) * 1000])

        ss_mleo_space_station.azimuth = param.azimuth
        ss_mleo_space_station.elevation = param.elevation

        ss_mleo_space_station.active = np.array([True])
        ss_mleo_space_station.tx_power = np.array(
            [param.tx_power_density + 10 * math.log10(param.bandwidth * 1e6) + 30])
        ss_mleo_space_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            ss_mleo_space_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R S.672":
            ss_mleo_space_station.antenna = np.array([AntennaS672(param)])
        elif param.antenna_pattern == "ITU-R S.1528 LEO":
            ss_mleo_space_station.antenna = np.array([AntennaS1528_LEO(param)])
        else:
            sys.stderr.write("ERROR\nInvalid SS-MLEO antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        ss_mleo_space_station.bandwidth = param.bandwidth
        ss_mleo_space_station.noise_temperature = param.noise_temperature
        ss_mleo_space_station.thermal_noise = -500
        ss_mleo_space_station.total_interference = -500

        return ss_mleo_space_station

    @staticmethod
    def generate_fss_space_station(param: ParametersFssSs):
        fss_space_station = StationManager(1)
        fss_space_station.station_type = StationType.FSS_SS

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A
        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km = (param.EARTH_RADIUS + param.altitude) / 1000
        dist_imt_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude) / 1000

        # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
        sat_lat_rad = param.lat_deg * np.pi / 180.
        imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
        x1 = dist_sat_centre_earth_km * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth_km * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
        z1 = dist_sat_centre_earth_km * np.sin(sat_lat_rad)

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = param.imt_lat_deg * np.pi / 180.
        fss_space_station.x = np.array([x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)]) * 1000
        fss_space_station.y = np.array([y1]) * 1000
        fss_space_station.height = np.array([(z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad)
                                              - dist_imt_centre_earth_km) * 1000])

        fss_space_station.azimuth = param.azimuth
        fss_space_station.elevation = param.elevation

        fss_space_station.active = np.array([True])
        fss_space_station.tx_power = np.array([param.tx_power_density + 10 * math.log10(param.bandwidth * 1e6) + 30])
        fss_space_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            fss_space_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R S.672":
            fss_space_station.antenna = np.array([AntennaS672(param)])
        elif param.antenna_pattern == "ITU-R S.1528":
            fss_space_station.antenna = np.array([AntennaS1528(param)])
        elif param.antenna_pattern == "FSS_SS":
            fss_space_station.antenna = np.array([AntennaFssSs(param)])
        else:
            sys.stderr.write("ERROR\nInvalid FSS SS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        fss_space_station.bandwidth = param.bandwidth
        fss_space_station.noise_temperature = param.noise_temperature
        fss_space_station.thermal_noise = -500
        fss_space_station.total_interference = -500

        return fss_space_station

    @staticmethod
    def generate_fss_earth_station(param: ParametersFssEs, random_number_gen: np.random.RandomState, *args):
        """
        Generates FSS Earth Station.

        Arguments:
            param: ParametersFssEs
            random_number_gen: np.random.RandomState
            topology (optional): Topology
        """
        if len(args):
            topology = args[0]

        fss_earth_station = StationManager(1)
        fss_earth_station.station_type = StationType.FSS_ES

        if param.location.upper() == "FIXED":
            fss_earth_station.x = np.array([param.x])
            fss_earth_station.y = np.array([param.y])
        elif param.location.upper() == "CELL":
            x, y, dummy1, dummy2 = StationFactory.get_random_position(1, topology, random_number_gen,
                                                                      param.min_dist_to_bs, True)
            fss_earth_station.x = np.array(x)
            fss_earth_station.y = np.array(y)
        elif param.location.upper() == "NETWORK":
            x, y, dummy1, dummy2 = StationFactory.get_random_position(1, topology, random_number_gen,
                                                                      param.min_dist_to_bs, False)
            fss_earth_station.x = np.array(x)
            fss_earth_station.y = np.array(y)
        elif param.location.upper() == "UNIFORM_DIST":
            # FSS ES is randomly (uniform) created inside a circle of radius
            # equal to param.max_dist_to_bs
            if param.min_dist_to_bs < 0:
                sys.stderr.write("ERROR\nInvalid minimum distance from FSS ES to BS: {}".format(param.min_dist_to_bs))
                sys.exit(1)
            while True:
                dist_x = random_number_gen.uniform(-param.max_dist_to_bs, param.max_dist_to_bs)
                dist_y = random_number_gen.uniform(-param.max_dist_to_bs, param.max_dist_to_bs)
                radius = np.sqrt(dist_x ** 2 + dist_y ** 2)
                if (radius > param.min_dist_to_bs) & (radius < param.max_dist_to_bs):
                    break
            fss_earth_station.x[0] = dist_x
            fss_earth_station.y[0] = dist_y
        else:
            sys.stderr.write("ERROR\nFSS-ES location type {} not supported".format(param.location))
            sys.exit(1)

        fss_earth_station.height = np.array([param.height])

        if param.azimuth.upper() == "RANDOM":
            fss_earth_station.azimuth = random_number_gen.uniform(-180., 180.)
        else:
            fss_earth_station.azimuth = float(param.azimuth)

        elevation = random_number_gen.uniform(param.elevation_min, param.elevation_max)
        fss_earth_station.elevation = np.array([elevation])

        fss_earth_station.active = np.array([True])
        fss_earth_station.tx_power = np.array([param.tx_power_density + 10 * math.log10(param.bandwidth * 1e6) + 30])
        fss_earth_station.rx_interference = -500

        if param.antenna_pattern.upper() == "OMNI":
            fss_earth_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern.upper() == "ITU-R S.1855":
            fss_earth_station.antenna = np.array([AntennaS1855(param)])
        elif param.antenna_pattern.upper() == "ITU-R S.465":
            fss_earth_station.antenna = np.array([AntennaS465(param)])
        elif param.antenna_pattern.upper() == "MODIFIED ITU-R S.465":
            fss_earth_station.antenna = np.array([AntennaModifiedS465(param)])
        elif param.antenna_pattern.upper() == "ITU-R S.580":
            fss_earth_station.antenna = np.array([AntennaS580(param)])
        else:
            sys.stderr.write("ERROR\nInvalid FSS ES antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        fss_earth_station.noise_temperature = param.noise_temperature
        fss_earth_station.bandwidth = np.array([param.bandwidth])
        fss_earth_station.noise_temperature = param.noise_temperature
        fss_earth_station.thermal_noise = -500
        fss_earth_station.total_interference = -500

        return fss_earth_station

    @staticmethod
    def generate_fs_station(param: ParametersFs):
        fs_station = StationManager(1)
        fs_station.station_type = StationType.FS

        if (param.channel_model == "P619"):
            # Coordinates according to ITU-R P619-1, Attachment A
            # calculate distances to the centre of the Earth
            dist_hibs_centre_earth_km = (param.EARTH_RADIUS + param.altitude) / 1000
            dist_system_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude) / 1000

            # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
            sat_lat_rad = param.hibs_lat_deg * np.pi / 180.
            imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
            x1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
            y1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
            z1 = dist_hibs_centre_earth_km * np.sin(sat_lat_rad)

            # rotate axis and calculate coordinates with origin at System
            sys_lat_rad = param.imt_lat_deg * np.pi / 180.
            fs_station.x = np.array([x1 * np.sin(sys_lat_rad) - z1 * np.cos(sys_lat_rad)]) * 1000
            fs_station.y = np.array([y1]) * 1000
            z2 = np.array(
                [(z1 * np.sin(sys_lat_rad) + x1 * np.cos(sys_lat_rad) - dist_system_centre_earth_km)]) * 1000
            fs_station.height = param.altitude - z2
        else:
            fs_station.x = np.array([param.x])
            fs_station.y = np.array([param.y])
            fs_station.height = np.array([param.height])
        if (param.distribution_enable == "ON"):
            if (param.distribution_type == "UNIFORM"):
                if (type(param.azimuth_distribution)) != list:
                    aux_azimuth = param.azimuth_distribution.split(',')
                    param.azimuth_distribution = [float(i) for i in aux_azimuth]
                    aux_elevation = param.elevation_distribution.split(',')
                    param.elevation_distribution = [float(i) for i in aux_elevation]
                param.azimuth = np.random.uniform(param.azimuth_distribution[0], param.azimuth_distribution[1])
                param.elevation = np.random.uniform(param.elevation_distribution[0], param.elevation_distribution[1])
                fs_station.azimuth = np.array([param.azimuth])
                fs_station.elevation = np.array([param.elevation])
            elif (param.distribution_type == "UNIFORM_NORMAL"):
                if (type(param.azimuth_distribution)) != list:
                    aux_azimuth = param.azimuth_distribution.split(',')
                    param.azimuth_distribution = [float(i) for i in aux_azimuth]
                    aux_elevation = param.elevation_distribution.split(',')
                    param.elevation_distribution = [float(i) for i in aux_elevation]
                param.azimuth = np.random.uniform(param.azimuth_distribution[0], param.azimuth_distribution[1])
                param.elevation = np.random.normal(param.elevation_distribution[0], param.elevation_distribution[1])
                fs_station.azimuth = np.array([param.azimuth])
                fs_station.elevation = np.array([param.elevation])
        else:
            fs_station.azimuth = np.array([param.azimuth])
            fs_station.elevation = np.array([param.elevation])

        fs_station.active = np.array([True])
        fs_station.tx_power = np.array([param.tx_power_density + 10 * math.log10(param.bandwidth * 1e6) + 30])
        fs_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            fs_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R F.699":
            fs_station.antenna = np.array([AntennaF699(param)])
        elif param.antenna_pattern == "ITU-R F.1336":
            fs_station.antenna = np.array([AntennaOmniF1336(param)])
        elif param.antenna_pattern == "BESSEL":
            fs_station.antenna = np.array([AntennaBessel(param)])
        elif param.antenna_pattern == "ITU-R F.1245":
            fs_station.antenna = np.array([AntennaF1245(param)])
        else:
            sys.stderr.write("ERROR\nInvalid FS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        fs_station.noise_temperature = param.noise_temperature
        fs_station.bandwidth = np.array([param.bandwidth])

        return fs_station

    @staticmethod
    def generate_haps(param: ParametersHaps, intersite_distance: int, random_number_gen: np.random.RandomState()):
        num_haps = 1
        haps = StationManager(num_haps)
        haps.station_type = StationType.HAPS

        #        d = intersite_distance
        #        h = (d/3)*math.sqrt(3)/2
        #        haps.x = np.array([0, 7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
        #        haps.y = np.array([0, 9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
        haps.x = np.array([0])
        haps.y = np.array([0])

        haps.height = param.altitude * np.ones(num_haps)

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
            sys.stderr.write("ERROR\nInvalid HAPS (airbone) antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        haps.bandwidth = np.array([param.bandwidth])

        return haps

    @staticmethod
    def generate_rns(param: ParametersRns, random_number_gen: np.random.RandomState()):
        num_rns = 1
        rns = StationManager(num_rns)
        rns.station_type = StationType.RNS

        rns.x = np.array([param.x])
        rns.y = np.array([param.y])
        rns.height = np.array([param.altitude])

        # minimum and maximum values for azimuth and elevation
        azimuth = np.array([-30, 30])
        elevation = np.array([-30, 5])

        rns.azimuth = 90 + (azimuth[1] - azimuth[0]) * random_number_gen.random_sample(num_rns) + azimuth[0]
        rns.elevation = (elevation[1] - elevation[0]) * random_number_gen.random_sample(num_rns) + elevation[0]

        rns.active = np.ones(num_rns, dtype=bool)

        if param.antenna_pattern == "OMNI":
            rns.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R M.1466":
            rns.antenna = np.array([AntennaM1466(param.antenna_gain, rns.azimuth, rns.elevation)])
        else:
            sys.stderr.write("ERROR\nInvalid RNS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        rns.bandwidth = np.array([param.bandwidth])
        rns.noise_temperature = param.noise_temperature
        rns.thermal_noise = -500
        rns.total_interference = -500
        rns.rx_interference = -500

        return rns

    @staticmethod
    def generate_arns(param: ParametersArns, random_number_gen: np.random.RandomState()):

        num_arns = 1
        arns = StationManager(num_arns)
        arns.station_type = StationType.ARNS

        if (param.channel_model == "P619"):
            # Coordinates according to  ITU-R P619-1, Attachment A
            # calculate distances to the centre of the Earth
            dist_hibs_centre_earth_km = (param.EARTH_RADIUS + param.altitude) / 1000
            dist_system_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude) / 1000

            # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
            sat_lat_rad = param.hibs_lat_deg * np.pi / 180.
            imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
            x1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
            y1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
            z1 = dist_hibs_centre_earth_km * np.sin(sat_lat_rad)

            # rotate axis and calculate coordinates with origin at System
            sys_lat_rad = param.imt_lat_deg * np.pi / 180.
            arns.x = np.array([x1 * np.sin(sys_lat_rad) - z1 * np.cos(sys_lat_rad)]) * 1000
            arns.y = np.array([y1]) * 1000
            z2 = np.array([(z1 * np.sin(sys_lat_rad) + x1 * np.cos(sys_lat_rad) - dist_system_centre_earth_km)]) * 1000
            arns.height = param.altitude - z2
            print(arns.height)
        else:
            arns.x = np.array([param.x])
            arns.y = np.array([param.y])
            arns.height = np.array([param.height])

        if (param.antenna_pattern=="PHASED ARRAY"):
            positions = [-180, -90, 0, 90]
            arns.azimuth = np.array([np.random.choice(positions)])
            arns.elevation = np.array([param.elevation])
            param.beamsteeringangle_az = np.random.uniform(-45,45)
            #param.beamsteeringangle_el = np.random.normal(-20,85)
            param.beamsteeringangle_el = np.random.uniform(-2.5,2.5)
            print(param.beamsteeringangle_az)
            print(param.beamsteeringangle_el)

        else:
            if (param.distribution_enable == "ON"):
                if (param.distribution_type == "UNIFORM"):
                    if (type(param.azimuth_distribution)) != list:
                        aux_azimuth = param.azimuth_distribution.split(',')
                        param.azimuth_distribution = [float(i) for i in aux_azimuth]
                        aux_elevation = param.elevation_distribution.split(',')
                        param.elevation_distribution = [float(i) for i in aux_elevation]
                    param.azimuth = np.random.uniform(param.azimuth_distribution[0], param.azimuth_distribution[1])
                    param.elevation = np.random.uniform(param.elevation_distribution[0], param.elevation_distribution[1])
                    arns.azimuth = np.array([param.azimuth])
                    arns.elevation = np.array([param.elevation])
                elif (param.distribution_type == "UNIFORM_NORMAL"):
                    if (type(param.azimuth_distribution)) != list:
                        aux_azimuth = param.azimuth_distribution.split(',')
                        param.azimuth_distribution = [float(i) for i in aux_azimuth]
                        aux_elevation = param.elevation_distribution.split(',')
                        param.elevation_distribution = [float(i) for i in aux_elevation]
                    param.azimuth = np.random.uniform(param.azimuth_distribution[0], param.azimuth_distribution[1])
                    param.elevation = np.random.normal(param.elevation_distribution[0], param.elevation_distribution[1])
                    arns.azimuth = np.array([param.azimuth])
                    arns.elevation = np.array([param.elevation])
            else:
                arns.azimuth = np.array([param.azimuth])
                arns.elevation = np.array([param.elevation])

        print(arns.azimuth)
        print(arns.elevation)

        arns.active = np.ones(num_arns, dtype=bool)

        if param.antenna_pattern == "OMNI":
            arns.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "COSSECANT SQUARED":
            arns.antenna = np.array([AntennaCossecantSquared(param)])
        elif param.antenna_pattern == "UNIFORM":
            arns.antenna = np.array([AntennaMeteorologicalRadarUniform(param)])
        elif param.antenna_pattern == "COSINE":
            arns.antenna = np.array([AntennaMeteorologicalRadarCosine1(param)])
        elif param.antenna_pattern == "PHASED ARRAY":
            arns.antenna = np.array([AntennaRadarPhasedArray(param)])
        else:
            sys.stderr.write("ERROR\nInvalid RNS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        arns.bandwidth = np.array([param.bandwidth])
        arns.noise_temperature = param.noise_temperature
        arns.thermal_noise = -500
        arns.total_interference = -500
        arns.rx_interference = -500
        return arns

    @staticmethod
    def generate_ras_station(param: ParametersRas):
        ras_station = StationManager(1)
        ras_station.station_type = StationType.RAS

        if (param.channel_model == "P619"):
            # Coordinates according to  ITU-R P619-1, Attachment A
            # calculate distances to the centre of the Earth
            dist_hibs_centre_earth_km = (param.EARTH_RADIUS + param.altitude) / 1000
            dist_system_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude) / 1000

            # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
            sat_lat_rad = param.hibs_lat_deg * np.pi / 180.
            imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
            x1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
            y1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
            z1 = dist_hibs_centre_earth_km * np.sin(sat_lat_rad)

            # rotate axis and calculate coordinates with origin at System
            sys_lat_rad = param.imt_lat_deg * np.pi / 180.
            ras_station.x = np.array([x1 * np.sin(sys_lat_rad) - z1 * np.cos(sys_lat_rad)]) * 1000
            ras_station.y = np.array([y1]) * 1000
            z2 = np.array([(z1 * np.sin(sys_lat_rad) + x1 * np.cos(sys_lat_rad) - dist_system_centre_earth_km)]) * 1000
            ras_station.height = param.altitude - z2
        else:
            ras_station.x = np.array([param.x])
            ras_station.y = np.array([param.y])
            ras_station.height = np.array([param.height])

        if (param.distribution_enable == "ON"):
            if (param.distribution_type == "UNIFORM"):
                if (type(param.azimuth_distribution)) != list:
                    aux_azimuth = param.azimuth_distribution.split(',')
                    param.azimuth_distribution = [int(i) for i in aux_azimuth]
                    aux_elevation = param.elevation_distribution.split(',')
                    param.elevation_distribution = [int(i) for i in aux_elevation]
                param.azimuth = np.random.uniform(param.azimuth_distribution[0], param.azimuth_distribution[1])
                param.elevation = np.random.uniform(param.elevation_distribution[0], param.elevation_distribution[1])
                ras_station.azimuth = np.array([param.azimuth])
                ras_station.elevation = np.array([param.elevation])
        else:
            ras_station.azimuth = np.array([param.azimuth])
            ras_station.elevation = np.array([param.elevation])

        ras_station.active = np.array([True])
        ras_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            ras_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
            ras_station.antenna[0].effective_area = param.SPEED_OF_LIGHT ** 2 / (
                4 * np.pi * (param.frequency * 1e6) ** 2)
        elif param.antenna_pattern == "ITU-R SA.509":
            ras_station.antenna = np.array([AntennaSA509(param)])
        else:
            sys.stderr.write("ERROR\nInvalid RAS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        ras_station.noise_temperature = np.array(param.antenna_noise_temperature + param.receiver_noise_temperature)
        ras_station.bandwidth = np.array(param.bandwidth)

        return ras_station

    @staticmethod
    def generate_eess_passive_sensor(param: ParametersEessPassive):
        eess_passive_sensor = StationManager(1)
        eess_passive_sensor.station_type = StationType.EESS_PASSIVE

        # incidence angle according to Rec. ITU-R RS.1861-0
        incidence_angle = math.degrees(math.asin(
            math.sin(math.radians(param.nadir_angle)) * (1 + (param.altitude / param.EARTH_RADIUS))))

        # distance to field of view centre according to Rec. ITU-R RS.1861-0
        distance = param.EARTH_RADIUS * \
                   math.sin(math.radians(incidence_angle - param.nadir_angle)) / \
                   math.sin(math.radians(param.nadir_angle))

        # Elevation at ground (centre of the footprint)
        theta_grd_elev = 90 - incidence_angle

        eess_passive_sensor.x = np.array([0])
        eess_passive_sensor.y = np.array([distance * math.cos(math.radians(theta_grd_elev))])
        eess_passive_sensor.height = np.array([distance * math.sin(math.radians(theta_grd_elev))])

        # Elevation and azimuth at sensor wrt centre of the footprint
        # It is assumed the sensor is at y-axis, hence azimuth is 270 deg
        eess_passive_sensor.azimuth = 270
        eess_passive_sensor.elevation = -theta_grd_elev

        eess_passive_sensor.active = np.array([True])
        eess_passive_sensor.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            eess_passive_sensor.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R RS.1813":
            eess_passive_sensor.antenna = np.array([AntennaRS1813(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9a":
            eess_passive_sensor.antenna = np.array([AntennaRS1861_9A(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9b":
            eess_passive_sensor.antenna = np.array([AntennaRS1861_9B(param)])
        elif param.antenna_pattern == "ITU-R RS.1861 9c":
            eess_passive_sensor.antenna = np.array([AntennaRS1861_9C()])
        else:
            sys.stderr.write("ERROR\nInvalid EESS PASSIVE antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        eess_passive_sensor.bandwidth = param.bandwidth
        # Noise temperature is not an input parameter for EESS passive.
        # It is included here to calculate the useless I/N values
        eess_passive_sensor.noise_temperature = 250
        eess_passive_sensor.thermal_noise = -500
        eess_passive_sensor.total_interference = -500

        return eess_passive_sensor

    @staticmethod
    def get_random_position(num_stas: int, topology: Topology,
                            random_number_gen: np.random.RandomState,
                            min_dist_to_bs=0, central_cell=False,
                            deterministic_cell=False, hibs=""):
        hexagon_radius = topology.intersite_distance / 3

        min_dist_ok = False

        while not min_dist_ok:
            # generate UE uniformly in a triangle
            x = random_number_gen.uniform(0, hexagon_radius * np.cos(np.pi / 6), num_stas)
            y = random_number_gen.uniform(0, hexagon_radius / 2, num_stas)

            invert_index = np.arctan(y / x) > np.pi / 6
            y[invert_index] = -(hexagon_radius / 2 - y[invert_index])
            x[invert_index] = (hexagon_radius * np.cos(np.pi / 6) - x[invert_index])

            if any(np.sqrt(x ** 2 + y ** 2) < min_dist_to_bs):
                min_dist_ok = False
            else:
                min_dist_ok = True

        # randomly choose an hextant
        hextant = random_number_gen.random_integers(0, 5, num_stas)
        hextant_angle = np.pi / 6 + np.pi / 3 * hextant

        old_x = x
        x = x * np.cos(hextant_angle) - y * np.sin(hextant_angle)
        y = old_x * np.sin(hextant_angle) + y * np.cos(hextant_angle)

        # randomly choose a cell
        if central_cell:
            central_cell_indices = np.where((topology.x == 0) & (topology.y == 0))
            cell = central_cell_indices[0][random_number_gen.random_integers(0, len(central_cell_indices[0]) - 1,
                                                                             num_stas)]
        elif deterministic_cell:
            num_bs = topology.num_base_stations
            stas_per_cell = num_stas / num_bs
            cell = np.repeat(np.arange(num_bs, dtype=int), stas_per_cell)
        else:
            num_bs = topology.num_base_stations
            cell = random_number_gen.random_integers(0, num_bs - 1, num_stas)

        cell_x = topology.x[cell]
        cell_y = topology.y[cell]

        x = x + cell_x + hexagon_radius * np.cos(topology.azimuth[cell] * np.pi / 180)
        y = y + cell_y + hexagon_radius * np.sin(topology.azimuth[cell] * np.pi / 180)

        x = list(x)
        y = list(y)

        # calculate UE azimuth wrt serving BS
        theta = np.arctan2(y - cell_y, x - cell_x)

        # calculate elevation angle
        # psi is the vertical angle of the UE wrt the serving BS
        distance = np.sqrt((cell_x - x) ** 2 + (cell_y - y) ** 2)

        return x, y, theta, distance

    @staticmethod
    def generate_imt_base_stations(param: ParametersImt,
                                   param_h: ParametersHibs,
                                   param_ant: ParametersAntennaImt,
                                   topology: Topology,
                                   random_number_gen: np.random.RandomState):
        # par = param_ant.get_antenna_parameters(StationType.IMT_BS)
        par = param_ant.get_antenna_parameters(StationType.IMT_BS)
        num_bs = topology.num_base_stations
        imt_base_stations = StationManager(num_bs)
        imt_base_stations.station_type = StationType.IMT_BS
        # now we set the coordinates
        imt_base_stations.x = topology.x
        imt_base_stations.y = topology.y
        imt_base_stations.azimuth = topology.azimuth
        bs_back_off_power = param_h.bs_backoff_power
        if param.topology == 'HIBS':
            imt_base_stations.elevation = topology.elevation
        else:
            imt_base_stations.elevation = -par.downtilt*np.ones(num_bs)

        if param.topology == 'INDOOR':
            imt_base_stations.height = topology.height
        elif param.topology == 'HIBS':
            imt_base_stations.height = topology.height * np.ones(num_bs)
        else:
            imt_base_stations.height = param.bs_height * np.ones(num_bs)

        imt_base_stations.active = random_number_gen.rand(num_bs) < param.bs_load_probability

        """ HIBs Base Station Conducted Power
                   ################################################################################
                   3 Cell Simulation
                   No modifications = bs_conducted_power
                   ################################################################################
                   7 Cell Simulation

                   1st Layer conducted power: Cell 0 has to power modification = bs_conducted_power
                                              Cell 1 to 6 you can add back-off power (Default = 3dB)
                   ################################################################################
                   19 Cell Simulation
                   No modifications = bs_conducted_power
                   ################################################################################
                   49 Cell Simulation (HIBS Cluster) - in development

                   1st Layer conducted power: Cell 0, 14, 21, 28, 35, 42 has to power modification = bs_conducted_power
                                              Other Cells - you can add back-off power (Default = 3dB)
        """

        if param.topology == 'HIBS':
            if num_bs == 1:
                imt_base_stations.tx_power = param_h.bs_conducted_power * np.ones(num_bs)

            elif num_bs == 3:
                imt_base_stations.tx_power = param_h.bs_conducted_power * np.ones(num_bs)

            elif num_bs == 7:
                imt_base_stations.tx_power =     [param_h.bs_conducted_power,
                                                  param_h.bs_conducted_power - bs_back_off_power,
                                                  param_h.bs_conducted_power - bs_back_off_power,
                                                  param_h.bs_conducted_power - bs_back_off_power,
                                                  param_h.bs_conducted_power - bs_back_off_power,
                                                  param_h.bs_conducted_power - bs_back_off_power,
                                                  param_h.bs_conducted_power - bs_back_off_power]

            elif num_bs == 19:
                imt_base_stations.tx_power = param_h.bs_conducted_power * np.ones(num_bs)

            elif num_bs == 49:
                imt_base_stations.tx_power = 7*[param_h.bs_conducted_power,
                                              param_h.bs_conducted_power - bs_back_off_power,
                                              param_h.bs_conducted_power - bs_back_off_power,
                                              param_h.bs_conducted_power - bs_back_off_power,
                                              param_h.bs_conducted_power - bs_back_off_power,
                                              param_h.bs_conducted_power - bs_back_off_power,
                                              param_h.bs_conducted_power - bs_back_off_power]
        else:
            imt_base_stations.tx_power = param.bs_conducted_power * np.ones(num_bs)

        #############################################################################

        imt_base_stations.rx_power = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.rx_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.ext_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.total_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])

        imt_base_stations.snr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.sinr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.sinr_ext = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.inr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])

        imt_base_stations.antenna = np.empty(num_bs, dtype=AntennaBeamformingImt)

        if param_ant.bs_antenna_type == 'BEAMFORMING':
            """ HIBs Base Station Antenna Array modification for 7 Cell Simulation
                            Configurate the antenna array number (row x colunn) per sector

                            ###################################################################################
                            For Cell 0 imt_base_stations.antenna[0]
                            change: par.bs_tx_n_rows = 2
                                    par.bs_tx_n_columns = 2
                            ###################################################################################
                            For other cells from 1 to 6
                            change:
                                    par.bs_tx_n_rows = 4
                                    par.bs_tx_n_columns = 2
                            ###################################################################################
             """

        if param.topology == 'HIBS':
            if num_bs == 7:
                par = param_ant
                par.bs_n_rows = 2
                par.bs_n_columns = 2
                par = par.get_antenna_parameters(StationType.IMT_BS)
                imt_base_stations.antenna[0] = AntennaBeamformingImt(par, imt_base_stations.azimuth[0],
                                                                     imt_base_stations.elevation[0])
                # For Cells from 1 to 6
                par = param_ant
                par.bs_n_rows = 4
                par.bs_n_columns = 2
                par = par.get_antenna_parameters(StationType.IMT_BS)
                for i in range(1, num_bs):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
            """ HIBs Base Station Antenna Array modification for 49 Cell Simulation (Cluster - in development)
                                        Configurate the antenna array number (row x colunn) per sector
            """
            if num_bs == 49:
                # Cell 0
                par = param_ant
                par.bs_n_rows = 2
                par.bs_n_columns = 2
                par = par.get_antenna_parameters(StationType.IMT_BS)
                imt_base_stations.antenna[0] = AntennaBeamformingImt(par, imt_base_stations.azimuth[0],
                                                                     imt_base_stations.elevation[0])
                imt_base_stations.antenna[7] = AntennaBeamformingImt(par, imt_base_stations.azimuth[7],
                                                                     imt_base_stations.elevation[7])
                imt_base_stations.antenna[14] = AntennaBeamformingImt(par, imt_base_stations.azimuth[14],
                                                                     imt_base_stations.elevation[14])
                imt_base_stations.antenna[21] = AntennaBeamformingImt(par, imt_base_stations.azimuth[21],
                                                                      imt_base_stations.elevation[21])
                imt_base_stations.antenna[28] = AntennaBeamformingImt(par, imt_base_stations.azimuth[28],
                                                                     imt_base_stations.elevation[28])
                imt_base_stations.antenna[35] = AntennaBeamformingImt(par, imt_base_stations.azimuth[35],
                                                                      imt_base_stations.elevation[35])
                imt_base_stations.antenna[42] = AntennaBeamformingImt(par, imt_base_stations.azimuth[42],
                                                                      imt_base_stations.elevation[42])
                # Other cells
                par = param_ant
                par.bs_n_rows = 4
                par.bs_n_columns = 2
                par = par.get_antenna_parameters(StationType.IMT_BS)
                #for i in range(1, num_bs):
                for i in range(1, 7):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(8, 14):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(15, 21):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(22, 28):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(29, 35):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(36, 42):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])
                for i in range(43, 49):
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])

        else:
            for i in range(num_bs):
                imt_base_stations.antenna[i] = \
                    imt_base_stations.antenna[i] = AntennaBeamformingImt(par, imt_base_stations.azimuth[i],
                                                                         imt_base_stations.elevation[i])

        # imt_base_stations.antenna = [AntennaOmni(0) for bs in range(num_bs)]
        imt_base_stations.bandwidth = param.bandwidth * np.ones(num_bs)
        imt_base_stations.center_freq = param.frequency * np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs_noise_figure * np.ones(num_bs)
        imt_base_stations.thermal_noise = -500 * np.ones(num_bs)

        if param.spectral_mask == "IMT-2020":
            imt_base_stations.spectral_mask = SpectralMaskImt(StationType.IMT_BS,
                                                              param.frequency,
                                                              param.bandwidth,
                                                              param.spurious_emissions,
                                                              scenario=param.topology)
        elif param.spectral_mask == "3GPP E-UTRA":
            imt_base_stations.spectral_mask = SpectralMask3Gpp(StationType.IMT_BS,
                                                               param.frequency,
                                                               param.bandwidth,
                                                               param.spurious_emissions)

        if param.topology == 'MACROCELL' or param.topology == 'HOTSPOT':
            imt_base_stations.intersite_dist = param.intersite_distance
        elif param.topology == 'HIBS':
            imt_base_stations.intersite_dist = param_h.intersite_distance

        return imt_base_stations


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    factory = StationFactory()

    class ParamsHibs(object):
        def __init__(self):
            self.cell_radius = 100000
            self.intersite_distance = self.cell_radius * np.sqrt(3)
            self.num_clusters = 1
            self.num_sectors = 7
            self.bs_height = 20000
            self.azimuth3 = '90,210,330'
            self.elevation3 = '-90,-90,-90'
            self.azimuth7 = '0,0,60,120,180,240,300'
            #self.elevation7 = '-90,-23,-23,-23,-23,-23,-23'
            self.elevation7 = '-90,-33,-33,-33,-33,-33,-33'
            self.azimuth19 = '0,15,30,45,75,90,105,135,150,165,195,210,225,255,270,285,315,330,345'
            self.elevation19 = '-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30'
            self.bs_conducted_power = 37
            self.bs_backoff_power = 3


    class ParamsAux(object):
        def __init__(self):
            self.ue_distribution_type = "UNIFORM"
            self.bs_height = 20000
            self.ue_height = 1.5
            self.ue_indoor_percent = 0
            self.ue_k = 3
            self.ue_k_m = 1
            self.bandwidth = 20
            self.ue_tx_power_control = 'ON'
            self.ue_noise_figure = -5
            self.topology = "HIBS"
            self.minimum_separation_distance_bs_ue = 0
            self.frequency = 2680
            self.spectral_mask = "3GPP E-UTRA"
            self.spurious_emissions = -13
            self.bs_load_probability = 1
            self.bs_noise_figure = 5


    params = ParamsAux()
    ant_param = ParametersAntennaImt()
    param_hibs_teste = ParamsHibs()
    ParametersHibs.bs_back_off_power = 3

    topology = TopologyHIBS(param_hibs_teste.intersite_distance, param_hibs_teste.cell_radius,
                            param_hibs_teste.num_clusters, param_hibs_teste.num_sectors, param_hibs_teste.bs_height,
                            param_hibs_teste.azimuth3, param_hibs_teste.azimuth7, param_hibs_teste.azimuth19,
                            param_hibs_teste.elevation3, param_hibs_teste.elevation7, param_hibs_teste.elevation19)

    topology.calculate_coordinates()

    ant_param.adjacent_antenna_model = "SINGLE_ELEMENT"  #
    ant_param.ue_normalization = False  #
    ant_param.bs_normalization = False  #

    ant_param.normalization = False
    ant_param.ue_normalization_data = ""
    ant_param.bs_normalization_data = ""

    ant_param.bs_antenna_type = "BEAMFORMING"
    ant_param.bf_enable = "OFF"
    ant_param.bs_element_pattern = "M2101"
    #    ant_param.bs_element_pattern = "F1336"
    ant_param.bs_element_max_g = 8
    ant_param.bs_element_phi_3db = 65
    ant_param.bs_element_theta_3db = 65
    ant_param.bs_element_am = 30
    ant_param.bs_element_sla_v = 30
    ant_param.bs_n_rows = 4
    ant_param.bs_n_columns = 4
    ant_param.bs_element_horiz_spacing = 0.5
    ant_param.bs_element_vert_spacing = 0.5
    ant_param.bs_multiplication_factor = 1  #
    ant_param.bs_minimum_array_gain = -200  #

    ant_param.bs_downtilt = 0

    ant_param.bs_rx_element_max_g = 8  #
    ant_param.bs_rx_element_phi_deg_3db = 65  #
    ant_param.bs_rx_element_theta_deg_3db = 65  #
    ant_param.bs_rx_element_am = 30  #
    ant_param.bs_rx_element_sla_v = 30  #
    ant_param.bs_rx_n_rows = 4  #
    ant_param.bs_rx_n_columns = 4  #
    ant_param.bs_rx_element_horiz_spacing = 0.5
    ant_param.bs_rx_element_vert_spacing = 0.5
    # ant_param.bs_antenna_type = 'M2101'

    ant_param.ue_antenna_type = 'BEAMFORMING'
    ant_param.ue_element_pattern = "FIXED"
    ant_param.ue_element_max_g = -3
    ant_param.ue_element_phi_3db = 90
    ant_param.ue_element_theta_3db = 90
    ant_param.ue_element_am = 23
    ant_param.ue_element_sla_v = 25
    ant_param.ue_n_rows = 4
    ant_param.ue_n_columns = 4
    ant_param.ue_element_horiz_spacing = 0.5
    ant_param.ue_element_vert_spacing = 0.5
    ant_param.ue_multiplication_factor = 1
    ant_param.ue_minimum_array_gain = -200

    imt_bs = factory.generate_imt_base_stations(params, param_hibs_teste, ant_param, topology, np.random.RandomState())
    imt_ue = factory.generate_imt_ue(params, ant_param, topology, np.random.RandomState())

    fig = plt.figure(figsize=(8, 8), linewidth=2, facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    topology.plot(ax)

    plt.axis('image')
    plt.title("HIBS Topology")
    plt.xlabel("x-coordinate [km]")
    plt.ylabel("y-coordinate [km]")

    plt.plot(imt_ue.x, imt_ue.y, "o", label="HIBS User Equipments", color='red')
    plt.legend()

    plt.tight_layout()
    plt.show()

    # from matplotlib import pyplot as plt
    #
    # # plot uniform distribution in macrocell scenario
    #
    # factory = StationFactory()
    # topology = TopologyMacrocell(1000, 1)
    # topology.calculate_coordinates()
    #
    # class ParamsAux(object):
    #     def __init__(self):
    #         self.spectral_mask = 'IMT-2020'
    #         self.frequency = 10000
    #         self.topology = 'MACROCELL'
    #         self.ue_distribution_type = "UNIFORM"
    #         self.bs_height = 30
    #         self.ue_height = 3
    #         self.ue_indoor_percent = 0
    #         self.ue_k = 3
    #         self.ue_k_m = 1
    #         self.bandwidth  = np.random.rand()
    #         self.ue_noise_figure = np.random.rand()
    #         self.minimum_separation_distance_bs_ue = 10
    #         self.spurious_emissions = -30
    #         self.intersite_distance = 1000
    #
    # params = ParamsAux()
    #
    # ant_param = ParametersAntennaImt()
    #
    # ant_param.adjacent_antenna_model = "SINGLE_ELEMENT"
    # ant_param.bs_element_pattern = "F1336"
    # ant_param.bs_element_max_g = 5
    # ant_param.bs_element_phi_3db = 65
    # ant_param.bs_element_theta_3db = 65
    # ant_param.bs_element_am = 30
    # ant_param.bs_element_sla_v = 30
    # ant_param.bs_n_rows = 8
    # ant_param.bs_n_columns = 8
    # ant_param.bs_element_horiz_spacing = 0.5
    # ant_param.bs_element_vert_spacing = 0.5
    # ant_param.bs_downtilt_deg = 10
    # ant_param.bs_multiplication_factor = 12
    # ant_param.bs_minimum_array_gain = -200
    #
    # ant_param.ue_element_pattern = "FIXED"
    # ant_param.ue_element_max_g = 5
    # ant_param.ue_element_phi_3db = 90
    # ant_param.ue_element_theta_3db = 90
    # ant_param.ue_element_am = 25
    # ant_param.ue_element_sla_v = 25
    # ant_param.ue_n_rows = 4
    # ant_param.ue_n_columns = 4
    # ant_param.ue_element_horiz_spacing = 0.5
    # ant_param.ue_element_vert_spacing = 0.5
    # ant_param.ue_multiplication_factor = 12
    # ant_param.ue_minimum_array_gain = -200
    #
    # ant_param.ue_normalization = False
    # ant_param.bs_normalization = False
    #
    # rnd = np.random.RandomState(1)
    #
    # imt_ue = factory.generate_imt_ue(params, ant_param, topology, rnd)
    #
    # fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')  # create a figure object
    # ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    #
    # topology.plot(ax)
    #
    # plt.axis('image')
    # plt.title("Macro cell topology")
    # plt.xlabel("x-coordinate [m]")
    # plt.ylabel("y-coordinate [m]")
    #
    # plt.plot(imt_ue.x, imt_ue.y, "r.")
    #
    # plt.tight_layout()
    # plt.show()
