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
from sharc.station_manager import StationManager
from sharc.spectral_mask_imt import SpectralMaskImt
from sharc.antenna.antenna import Antenna
from sharc.antenna.antenna_fss_ss import AntennaFssSs
from sharc.antenna.antenna_omni import AntennaOmni
from sharc.antenna.antenna_f699 import AntennaF699
from sharc.antenna.antenna_f1891 import AntennaF1891
from sharc.antenna.antenna_m1466 import AntennaM1466
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_s580 import AntennaS580
from sharc.antenna.antenna_s672 import AntennaS672
from sharc.antenna.antenna_s1528 import AntennaS1528
from sharc.antenna.antenna_s1855 import AntennaS1855
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.topology.topology import Topology

class StationFactory(object):

    @staticmethod
    def generate_imt_base_stations(param: ParametersImt,
                                   param_ant: ParametersAntennaImt,
                                   topology: Topology):
        num_bs = topology.num_base_stations
        imt_base_stations = StationManager(num_bs)
        imt_base_stations.station_type = StationType.IMT_BS
        # now we set the coordinates
        imt_base_stations.x = topology.x
        imt_base_stations.y = topology.y
        imt_base_stations.azimuth = topology.azimuth
        imt_base_stations.elevation = topology.elevation
        imt_base_stations.height = param.bs_height*np.ones(num_bs)
        imt_base_stations.indoor = topology.indoor
        imt_base_stations.active = np.random.rand(num_bs) < param.bs_load_probability
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
        
        imt_base_stations.spectral_mask = SpectralMaskImt(StationType.IMT_BS,param.frequency,\
                                                          param.bandwidth,scenario = param.topology)
        
        return imt_base_stations


    @staticmethod
    def generate_imt_ue(param: ParametersImt,
                        param_ant: ParametersAntennaImt,
                        topology: Topology) -> StationManager:
        
        if param.topology == "INDOOR":
            return StationFactory.generate_imt_ue_indoor(param, param_ant, topology)
        else:
            return StationFactory.generate_imt_ue_outdoor(param, param_ant, topology)

            
    @staticmethod
    def generate_imt_ue_outdoor(param: ParametersImt,
                                param_ant: ParametersAntennaImt,
                                topology: Topology) -> StationManager:           
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k*param.ue_k_m
        num_ue = num_bs*num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE
        ue_x = list()
        ue_y = list()

        # The Rayleigh and Normal distribution parameters (mean, scale and cutoff)
        # were agreed in TG 5/1 meeting (May 2017).

        if param.ue_distribution_distance == "RAYLEIGH":
            # For the distance between UE and BS, it is desired that 99% of UE's
            # are located inside the [soft] cell edge, i.e. Prob(d<d_edge) = 99%.
            # Since the distance is modeled by a random variable with Rayleigh
            # distribution, we use the quantile function to find that
            # sigma = distance/3.0345. So we always distibute UE's in order to meet
            # the requirement Prob(d<d_edge) = 99% for a given cell radius.
            radius_scale = topology.cell_radius/3.0345
            radius = np.random.rayleigh(radius_scale, num_ue)
        elif param.ue_distribution_distance == "UNIFORM":
            radius = topology.cell_radius*np.random.random(num_ue)
        else:
            sys.stderr.write("ERROR\nInvalid UE distance distribution: " + param.ue_distribution_distance)
            sys.exit(1)

        if param.ue_distribution_azimuth == "NORMAL":
            # In case of the angles, we generate N times the number of UE's because
            # the angle cutoff will discard 5% of the terminals whose angle is
            # outside the angular sector defined by [-60, 60]. So, N = 1.4 seems to
            # be a safe choice.
            N = 1.4
            angle_scale = 30
            angle_mean = 0
            angle_n = np.random.normal(angle_mean, angle_scale, int(N*num_ue))

            angle_cutoff = 60
            idx = np.where((angle_n < angle_cutoff) & (angle_n > -angle_cutoff))[0][:num_ue]
            angle = angle_n[idx]
        elif param.ue_distribution_azimuth == "UNIFORM":
            azimuth_range = (-60, 60)
            angle = (azimuth_range[1] - azimuth_range[0])*np.random.random(num_ue) + azimuth_range[0]
        else:
            sys.stderr.write("ERROR\nInvalid UE azimuth distribution: " + param.ue_distribution_distance)
            sys.exit(1)

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0])*np.random.random(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        #azimuth = np.zeros(num_ue)
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0])*np.random.random(num_ue) + elevation_range[0]

        for bs in range(num_bs):
            idx = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
            # theta is the horizontal angle of the UE wrt the serving BS
            theta = topology.azimuth[bs] + angle[idx]
            # calculate UE position in x-y coordinates
            x = topology.x[bs] + radius[idx]*np.cos(np.radians(theta))
            y = topology.y[bs] + radius[idx]*np.sin(np.radians(theta))
            ue_x.extend(x)
            ue_y.extend(y)
            # calculate UE azimuth wrt serving BS
            imt_ue.azimuth[idx] = (azimuth[idx] + theta + 180)%360

            # calculate elevation angle
            # psi is the vertical angle of the UE wrt the serving BS
            distance = np.sqrt((topology.x[bs] - x)**2 + (topology.y[bs] - y)**2)
            psi = np.degrees(np.arctan((param.bs_height - param.ue_height)/distance))
            imt_ue.elevation[idx] = elevation[idx] + psi

        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.height = param.ue_height*np.ones(num_ue)
        imt_ue.indoor = np.random.random(num_ue) <= (param.ue_indoor_percent/100)
        imt_ue.tx_power = param.ue_conducted_power*np.ones(num_ue)
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
        
        imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,param.frequency,\
                                               param.bandwidth,scenario = "OUTDOOR")
        
        return imt_ue

        
    @staticmethod
    def generate_imt_ue_indoor(param: ParametersImt,
                               param_ant: ParametersAntennaImt,
                               topology: Topology) -> StationManager:           
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k*param.ue_k_m
        num_ue = num_bs*num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE
        ue_x = list()
        ue_y = list()
        
        # initially set all UE's as indoor
        imt_ue.indoor = np.ones(num_ue, dtype=bool)

        # Calculate UE pointing
        azimuth_range = (-60, 60)
        azimuth = (azimuth_range[1] - azimuth_range[0])*np.random.random(num_ue) + azimuth_range[0]
        # Remove the randomness from azimuth and you will have a perfect pointing
        #azimuth = np.zeros(num_ue)
        elevation_range = (-90, 90)
        elevation = (elevation_range[1] - elevation_range[0])*np.random.random(num_ue) + elevation_range[0]
        
        delta_x = (topology.b_w/math.sqrt(topology.ue_indoor_percent) - topology.b_w)/2
        delta_y = (topology.b_d/math.sqrt(topology.ue_indoor_percent) - topology.b_d)/2

        for bs in range(num_bs):
            idx = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
            if bs % 3 == 0:
                x_min = topology.x[bs] - topology.cell_radius - delta_x
                x_max = topology.x[bs] + topology.cell_radius
            if bs % 3 == 1:
                x_min = topology.x[bs] - topology.cell_radius
                x_max = topology.x[bs] + topology.cell_radius
            if bs % 3 == 2:
                x_min = topology.x[bs] - topology.cell_radius
                x_max = topology.x[bs] + topology.cell_radius + delta_x
            y_min = topology.y[bs] - topology.b_d/2 - delta_y
            y_max = topology.y[bs] + topology.b_d/2 + delta_y
            x = (x_max - x_min)*np.random.random(num_ue_per_bs) + x_min
            y = (y_max - y_min)*np.random.random(num_ue_per_bs) + y_min
            ue_x.extend(x)
            ue_y.extend(y)
        
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
            if bs % 3 == 0:
                out = (x < topology.x[bs] - topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            if bs % 3 == 1:
                out = (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            if bs % 3 == 2:
                out = (x > topology.x[bs] + topology.cell_radius) | \
                      (y > topology.y[bs] + topology.b_d/2) | \
                      (y < topology.y[bs] - topology.b_d/2)
            imt_ue.indoor[idx] = ~ out
                
        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)

        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.height = param.ue_height*np.ones(num_ue)
        imt_ue.tx_power = param.ue_conducted_power*np.ones(num_ue)
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
        
        imt_ue.spectral_mask = SpectralMaskImt(StationType.IMT_UE,param.frequency,\
                                               param.bandwidth,scenario = "INDOOR")
        imt_ue.spectral_mask.set_mask()
        
        return imt_ue
        

    @staticmethod
    def generate_system(parameters: Parameters):
        if parameters.general.system == "FSS_ES":
            return StationFactory.generate_fss_earth_station(parameters.fss_es)
        elif parameters.general.system == "FSS_SS":
            return StationFactory.generate_fss_space_station(parameters.fss_ss)
        elif parameters.general.system == "FS":
            return StationFactory.generate_fs_station(parameters.fs)
        elif parameters.general.system == "HAPS":
            return StationFactory.generate_haps(parameters.haps, parameters.imt.intersite_distance)
        elif parameters.general.system == "RNS":
            return StationFactory.generate_rns(parameters.rns)
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
    def generate_fss_earth_station(param: ParametersFssEs):
        fss_earth_station = StationManager(1)
        fss_earth_station.station_type = StationType.FSS_ES

        fss_earth_station.x = np.array([param.x])
        fss_earth_station.y = np.array([param.y])
        fss_earth_station.height = np.array([param.height])

        fss_earth_station.azimuth = np.array([param.azimuth])
        fss_earth_station.elevation = np.array([param.elevation])

        fss_earth_station.active = np.array([True])
        fss_earth_station.tx_power = np.array([param.tx_power_density + 10*math.log10(param.bandwidth*1e6) + 30])
        fss_earth_station.rx_interference = -500

        if param.antenna_pattern == "OMNI":
            fss_earth_station.antenna = np.array([AntennaOmni(param.antenna_gain)])
        elif param.antenna_pattern == "ITU-R S.1855":
            fss_earth_station.antenna = np.array([AntennaS1855(param)])
        elif param.antenna_pattern == "ITU-R S.465":
            fss_earth_station.antenna = np.array([AntennaS465(param)])
        elif param.antenna_pattern == "ITU-R S.580":
            fss_earth_station.antenna = np.array([AntennaS580(param)])  
        else:
            sys.stderr.write("ERROR\nInvalid FSS ES antenna pattern: " + param.antenna_pattern)
            sys.exit(1)

        fss_earth_station.bandwidth = np.array([param.bandwidth])
        
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
            
        fs_station.bandwidth = np.array([param.bandwidth])
        
        return fs_station
        
        
    @staticmethod
    def generate_haps(param: ParametersHaps, intersite_distance: int):
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

        #haps.azimuth = param.azimuth
        #haps.elevation = param.elevation

        elev_max = 68.19 # corresponds to 50 km radius and 20 km altitude
        haps.azimuth = 360 * np.random.random(num_haps)
        haps.elevation = ((270 + elev_max) - (270 - elev_max)) * np.random.random(num_haps) + (270 - elev_max)
        
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
    def generate_rns(param: ParametersRns):
        num_rns = 1
        rns = StationManager(num_rns)
        rns.station_type = StationType.RNS

        rns.x = np.array([param.x])
        rns.y = np.array([param.y])
        rns.height = np.array([param.altitude])
        
        # minimum and maximum values for azimuth and elevation
        azimuth = np.array([-30, 30])
        elevation = np.array([-30, 5])

        rns.azimuth = 90 + (azimuth[1] - azimuth[0]) * np.random.random(num_rns) + azimuth[0]
        rns.elevation = (elevation[1] - elevation[0]) * np.random.random(num_rns) + elevation[0]
        
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
        
        