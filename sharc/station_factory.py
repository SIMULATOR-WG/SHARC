# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
"""

import numpy as np
import sys
import math

from sharc.support.enumerations import StationType
from sharc.parameters.parameters import Parameters
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fs import ParametersFs
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.parameters.parameters_haps import ParametersHaps
from sharc.parameters.parameters_rns import ParametersRns
from sharc.parameters.parameters_ras import ParametersRas
from sharc.station_manager import StationManager
from sharc.spectral_mask_imt import SpectralMaskImt
from sharc.antenna.antenna import Antenna
from sharc.antenna.antenna_fss_ss import AntennaFssSs
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.antenna.antenna_f699 import AntennaF699
from sharc.antenna.antenna_f1891 import AntennaF1891
from sharc.antenna.antenna_m1466 import AntennaM1466
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_modified_s465 import AntennaModifiedS465
from sharc.antenna.antenna_s580 import AntennaS580
from sharc.antenna.antenna_s672 import AntennaS672
from sharc.antenna.antenna_s1528 import AntennaS1528
from sharc.antenna.antenna_s1855 import AntennaS1855
from sharc.antenna.antenna_sa509 import AntennaSA509
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.spectral_mask_3gpp import SpectralMask3Gpp


class StationFactory(object):

    @staticmethod
    def generate_imt_base_stations(param: ParametersImt,
                                   param_ant: ParametersAntennaImt,
                                   topology: Topology,
                                   random_number_gen: np.random.RandomState):
        num_bs = topology.num_base_stations
        imt_base_stations = StationManager(num_bs)
        imt_base_stations.station_type = StationType.IMT_BS
        # now we set the coordinates
        imt_base_stations.x = topology.x
        imt_base_stations.y = topology.y
        imt_base_stations.azimuth = topology.azimuth
        imt_base_stations.elevation = topology.elevation
        if param.topology == 'INDOOR':
            imt_base_stations.height = topology.height
        else:
            imt_base_stations.height = param.bs_height*np.ones(num_bs)

        imt_base_stations.active = random_number_gen.rand(num_bs) < param.bs_load_probability
        imt_base_stations.tx_power = param.bs_conducted_power*np.ones(num_bs)
        imt_base_stations.rx_power = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.rx_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.ext_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.total_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])

        imt_base_stations.snr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.sinr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.sinr_ext = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.inr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])

        imt_base_stations.antenna = np.empty(num_bs, dtype=AntennaBeamformingImt)
        par = param_ant.get_antenna_parameters("BS", "RX")

        for i in range(num_bs):
            imt_base_stations.antenna[i] = \
            AntennaBeamformingImt(par, imt_base_stations.azimuth[i],\
                                  imt_base_stations.elevation[i])

        #imt_base_stations.antenna = [AntennaOmni(0) for bs in range(num_bs)]
        imt_base_stations.bandwidth = param.bandwidth*np.ones(num_bs)
        imt_base_stations.center_freq = param.frequency*np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs_noise_figure*np.ones(num_bs)
        imt_base_stations.thermal_noise = -500*np.ones(num_bs)

        if param.spectral_mask == "ITU 265-E":
            imt_base_stations.spectral_mask = SpectralMaskImt(StationType.IMT_BS,param.frequency,\
                                                              param.bandwidth,scenario = param.topology)
        elif param.spectral_mask == "3GPP 36.104":
            imt_base_stations.spectral_mask = SpectralMask3Gpp(StationType.IMT_BS,param.frequency,\
                                                               param.bandwidth)

        return imt_base_stations

    @staticmethod
    def generate_imt_ue(param: ParametersImt,
                        param_ant: ParametersAntennaImt,
                        topology: Topology,
                        random_number_gen: np.random.RandomState)-> StationManager:

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
        num_ue_per_bs = param.ue_k*param.ue_k_m

        num_ue = num_bs * num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE

        ue_x = list()
        ue_y = list()

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0])*random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0])*random_number_gen.random_sample(num_ue) + \
                    elevation_range[0]

        if param.ue_distribution_type.upper() == "UNIFORM":

            if not (type(topology) is TopologyMacrocell):
                sys.stderr.write("ERROR\nUniform UE distribution is currently supported only with Macrocell topology")
                sys.exit(1)

            [ue_x, ue_y, theta, distance] = StationFactory.get_random_position(num_ue, topology, random_number_gen,
                                                                               param.minimum_separation_distance_bs_ue )
            psi = np.degrees(np.arctan((param.bs_height - param.ue_height) / distance))

            imt_ue.azimuth = (azimuth + theta + np.pi/2)
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
        imt_ue.height = param.ue_height*np.ones(num_ue)
        imt_ue.indoor = random_number_gen.random_sample(num_ue) <= (param.ue_indoor_percent/100)
        imt_ue.rx_interference = -500*np.ones(num_ue)
        imt_ue.ext_interference = -500*np.ones(num_ue)

        # TODO: this piece of code works only for uplink
        par = param_ant.get_antenna_parameters("UE","TX")
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(par, imt_ue.azimuth[i],
                                                           imt_ue.elevation[i])

        #imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth*np.ones(num_ue)
        imt_ue.center_freq = param.frequency*np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure*np.ones(num_ue)

        if param.spectral_mask == "ITU 265-E":
            imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,param.frequency,\
                                                   param.bandwidth,scenario = "OUTDOOR")

        elif param.spectral_mask == "3GPP 36.104":
            imt_ue.spectral_mask = SpectralMask3Gpp(StationType.IMT_UE,param.frequency,\
                                                   param.bandwidth)

        imt_ue.spectral_mask.set_mask()

        return imt_ue


    @staticmethod
    def generate_imt_ue_indoor(param: ParametersImt,
                               param_ant: ParametersAntennaImt,
                               random_number_gen: np.random.RandomState,
                               topology: Topology) -> StationManager:
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k*param.ue_k_m
        num_ue = num_bs*num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE
        ue_x = list()
        ue_y = list()
        ue_z = list()

        # initially set all UE's as indoor
        imt_ue.indoor = np.ones(num_ue, dtype=bool)

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0])*random_number_gen.random_sample(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        #azimuth = np.zeros(num_ue)
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0])*random_number_gen.random_sample(num_ue) + elevation_range[0]

        delta_x = (topology.b_w/math.sqrt(topology.ue_indoor_percent) - topology.b_w)/2
        delta_y = (topology.b_d/math.sqrt(topology.ue_indoor_percent) - topology.b_d)/2

        for bs in range(num_bs):
            idx = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
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
                y_min = topology.y[bs] - topology.b_d/2 - delta_y
                y_max = topology.y[bs] + topology.b_d/2 + delta_y
            # Higher floors
            else:
                y_min = topology.y[bs] - topology.b_d/2
                y_max = topology.y[bs] + topology.b_d/2
                
            x = (x_max - x_min)*random_number_gen.random_sample(num_ue_per_bs) + x_min
            y = (y_max - y_min)*random_number_gen.random_sample(num_ue_per_bs) + y_min
            z = [topology.height[bs] - topology.b_h + param.ue_height for k in range(num_ue_per_bs)]
            ue_x.extend(x)
            ue_y.extend(y)
            ue_z.extend(z)

            # theta is the horizontal angle of the UE wrt the serving BS
            theta = np.degrees(np.arctan2(y - topology.y[bs], x - topology.x[bs]))
            # calculate UE azimuth wrt serving BS
            imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180)%360

            # calculate elevation angle
            # psi is the vertical angle of the UE wrt the serving BS
            distance = np.sqrt((topology.x[bs] - x)**2 + (topology.y[bs] - y)**2)
            psi = np.degrees(np.arctan((param.bs_height - param.ue_height)/distance))
            imt_ue.elevation[idx] = elevation[idx] + psi

            # check if UE is indoor
            if bs % topology.num_cells == 0:
                out = (x < topology.x[bs] - topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            elif bs % topology.num_cells == topology.num_cells - 1:
                out = (x > topology.x[bs] + topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            else:
                out = (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            imt_ue.indoor[idx] = ~ out

        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)
        imt_ue.height = np.array(ue_z)

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.rx_interference = -500*np.ones(num_ue)
        imt_ue.ext_interference = -500*np.ones(num_ue)

        # TODO: this piece of code works only for uplink
        par = param_ant.get_antenna_parameters("UE","TX")
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(par, imt_ue.azimuth[i],
                                                         imt_ue.elevation[i])

        #imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        imt_ue.bandwidth = param.bandwidth*np.ones(num_ue)
        imt_ue.center_freq = param.frequency*np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure*np.ones(num_ue)

        if param.spectral_mask == "ITU 265-E":
            imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,param.frequency,\
                                                   param.bandwidth,scenario = "INDOOR")

        elif param.spectral_mask == "3GPP 36.104":
            imt_ue.spectral_mask = SpectralMask3Gpp(StationType.IMT_UE,param.frequency,\
                                                   param.bandwidth)

        imt_ue.spectral_mask.set_mask()

        return imt_ue


    @staticmethod
    def generate_system(parameters: Parameters, topology: Topology, random_number_gen: np.random.RandomState ):
        if parameters.general.system == "FSS_ES":
            return StationFactory.generate_fss_earth_station(parameters.fss_es, random_number_gen, topology)
        elif parameters.general.system == "FSS_SS":
            return StationFactory.generate_fss_space_station(parameters.fss_ss)
        elif parameters.general.system == "FS":
            return StationFactory.generate_fs_station(parameters.fs)
        elif parameters.general.system == "HAPS":
            return StationFactory.generate_haps(parameters.haps, parameters.imt.intersite_distance, random_number_gen)
        elif parameters.general.system == "RNS":
            return StationFactory.generate_rns(parameters.rns, random_number_gen)
        elif parameters.general.system == "RAS":
            return StationFactory.generate_ras_station(parameters.ras)
        else:
            sys.stderr.write("ERROR\nInvalid system: " + parameters.general.system)
            sys.exit(1)


    @staticmethod
    def generate_fss_space_station(param: ParametersFssSs):
        fss_space_station = StationManager(1)
        fss_space_station.station_type = StationType.FSS_SS

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km = (param.EARTH_RADIUS + param.altitude)/1000
        dist_imt_centre_earth_km = (param.EARTH_RADIUS + param.imt_altitude)/1000

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
        fss_space_station.tx_power = np.array([param.tx_power_density + 10*math.log10(param.bandwidth*1e6) + 30])
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
        if len(args): topology = args[0]

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
            dist = random_number_gen.uniform( param.min_dist_to_bs, param.max_dist_to_bs)
            angle = random_number_gen.uniform(-np.pi, np.pi)
            fss_earth_station.x[0] = np.array(dist * np.cos(angle))
            fss_earth_station.y[0] = np.array(dist * np.sin(angle))
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
        fss_earth_station.tx_power = np.array([param.tx_power_density + 10*math.log10(param.bandwidth*1e6) + 30])
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

        fs_station.x = np.array([param.x])
        fs_station.y = np.array([param.y])
        fs_station.height = np.array([param.height])

        fs_station.azimuth = np.array([param.azimuth])
        fs_station.elevation = np.array([param.elevation])

        fs_station.active = np.array([True])
        fs_station.tx_power = np.array([param.tx_power_density + 10*math.log10(param.bandwidth*1e6) + 30])
        fs_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            fs_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R F.699":
            fs_station.antenna = np.array([AntennaF699(param)])
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

        elev_max = 68.19 # corresponds to 50 km radius and 20 km altitude
        haps.azimuth = 360 * random_number_gen.random_sample(num_haps)
        haps.elevation = ((270 + elev_max) - (270 - elev_max)) * random_number_gen.random_sample(num_haps) + \
                         (270 - elev_max)

        haps.active = np.ones(num_haps, dtype = bool)

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

        rns.active = np.ones(num_rns, dtype = bool)

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
    def generate_ras_station(param: ParametersRas):
        ras_station = StationManager(1)
        ras_station.station_type = StationType.RAS

        ras_station.x = np.array([param.x])
        ras_station.y = np.array([param.y])
        ras_station.height = np.array([param.height])

        ras_station.azimuth = np.array([param.azimuth])
        ras_station.elevation = np.array([param.elevation])

        ras_station.active = np.array([True])
        ras_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            ras_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R SA.509":
            ras_station.antenna = np.array([AntennaSA509(param)])
        else:
            sys.stderr.write("ERROR\nInvalid RAS antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        ras_station.noise_temperature = np.array(param.antenna_noise_temperature + \
                                                  param.receiver_noise_temperature)
        ras_station.bandwidth = np.array(param.bandwidth)

        return ras_station

    @staticmethod
    def get_random_position( num_stas: int, topology: Topology,
                             random_number_gen: np.random.RandomState,
                             min_dist_to_bs = 0, central_cell = False ):
        hexagon_radius = topology.intersite_distance / 3

        min_dist_ok = False

        while not min_dist_ok:
            # generate UE uniformly in a triangle
            x = random_number_gen.uniform(0, hexagon_radius * np.cos(np.pi / 6), num_stas)
            y = random_number_gen.uniform(0, hexagon_radius / 2, num_stas)

            invert_index = np.arctan(y / x) > np.pi / 6
            y[invert_index] = -(hexagon_radius / 2 - y[invert_index])
            x[invert_index] = (hexagon_radius * np.cos(np.pi / 6) - x[invert_index])

            if any (np.sqrt(x**2 + y**2) <  min_dist_to_bs):
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


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    # plot uniform distribution in macrocell scenario

    factory = StationFactory()
    topology = TopologyMacrocell(1000, 1)
    topology.calculate_coordinates()

    class ParamsAux(object):
        def __init__(self):
            self.ue_distribution_type = "UNIFORM"
            self.bs_height = 30
            self.ue_height = 3
            self.ue_indoor_percent = 0
            self.ue_k = 3
            self.ue_k_m = 20
            self.bandwidth  = np.random.rand()
            self.ue_noise_figure = np.random.rand()

    params = ParamsAux()

    ant_param = ParametersAntennaImt()

    ant_param.bs_element_pattern = "F1336"
    ant_param.bs_tx_element_max_g = 5
    ant_param.bs_tx_element_phi_deg_3db = 65
    ant_param.bs_tx_element_theta_deg_3db = 65
    ant_param.bs_tx_element_am = 30
    ant_param.bs_tx_element_sla_v = 30
    ant_param.bs_tx_n_rows = 8
    ant_param.bs_tx_n_columns = 8
    ant_param.bs_tx_element_horiz_spacing = 0.5
    ant_param.bs_tx_element_vert_spacing = 0.5
    ant_param.bs_downtilt_deg = 10

    ant_param.ue_element_pattern = "FIXED"
    ant_param.ue_tx_element_max_g = 5
    ant_param.ue_tx_element_phi_deg_3db = 90
    ant_param.ue_tx_element_theta_deg_3db = 90
    ant_param.ue_tx_element_am = 25
    ant_param.ue_tx_element_sla_v = 25
    ant_param.ue_tx_n_rows = 4
    ant_param.ue_tx_n_columns = 4
    ant_param.ue_tx_element_horiz_spacing = 0.5
    ant_param.ue_tx_element_vert_spacing = 0.5

    imt_ue = factory.generate_imt_ue(params, ant_param, topology)

    fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    topology.plot(ax)

    plt.axis('image')
    plt.title("Macro cell topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")

    plt.plot(imt_ue.x, imt_ue.y, "*")

    plt.tight_layout()
    plt.show()
