# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
"""

import numpy as np

from sharc.support.enumerations import StationType
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_manager import StationManager
from sharc.antenna.antenna_omni import AntennaOmni
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
        
        imt_base_stations.active = np.random.rand(num_bs) < param.bs_load_probability
        imt_base_stations.tx_power = param.bs_tx_power*np.ones(num_bs)
        imt_base_stations.rx_power = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.rx_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.total_interference = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        
        imt_base_stations.snr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        imt_base_stations.sinr = dict([(bs, -500 * np.ones(param.ue_k)) for bs in range(num_bs)])
        
        #imt_base_stations.antenna = [AntennaOmni(0) for bs in range(num_bs)]
        imt_base_stations.antenna = np.empty(num_bs, dtype=AntennaBeamformingImt)
        par = param_ant.get_antenna_parameters("BS", "RX")
        
        for i in range(num_bs):
            imt_base_stations.antenna[i] = \
            AntennaBeamformingImt(par, imt_base_stations.azimuth[i],\
                                  imt_base_stations.elevation[i])
        
        imt_base_stations.bandwidth = param.bandwidth*np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs_noise_figure*np.ones(num_bs)
        imt_base_stations.thermal_noise = -500*np.ones(num_bs)
        return imt_base_stations

        
    @staticmethod
    def generate_imt_ue(param: ParametersImt, 
                        param_ant: ParametersAntennaImt,
                        topology: Topology):
        num_bs = topology.num_base_stations
        num_ue_per_bs = param.ue_k*param.ue_k_m
        num_ue = num_bs*num_ue_per_bs

        imt_ue = StationManager(num_ue)
        imt_ue.station_type = StationType.IMT_UE
        ue_x = list()
        ue_y = list()
        
        # The Rayleigh and Normal distribution parameters (mean, scale and cutoff)
        # were agreed in TG 5/1 meeting (May 2017).
        
        # For the distance between UE and BS, it is desired that 99% of UE's 
        # are located inside the [soft] cell edge, i.e. Prob(d<d_edge) = 99%.
        # Since the distance is modeled by a random variable with Rayleigh
        # distribution, we use the quantile function to find that 
        # sigma = distance/3.0345. So we always distibute UE's in order to meet
        # the requirement Prob(d<d_edge) = 99% for a given cell radius.
        radius_scale = topology.cell_radius/3.0345
        radius = np.random.rayleigh(radius_scale, num_ue)
        
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
        imt_ue.tx_power = param.ue_tx_power*np.ones(num_ue)
        imt_ue.rx_interference = -500*np.ones(num_ue)

        #imt_ue.antenna = [AntennaOmni(0) for bs in range(num_ue)]
        # TODO: this piece of code works only for uplink
        par = param_ant.get_antenna_parameters("UE","TX")
        for i in range(num_ue):
            imt_ue.antenna[i] = AntennaBeamformingImt(par, imt_ue.azimuth[i], 
                                                         imt_ue.elevation[i])
            
        imt_ue.bandwidth = param.bandwidth*np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure*np.ones(num_ue)
        return imt_ue

        
    @staticmethod
    def generate_fss_stations(param: ParametersFss):
        satellite_stations = StationManager(1)
        satellite_stations.station_type = StationType.FSS_SS

        # now we set the coordinates according to
        # ITU-R P619-1, Attachment A

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth = param.EARTH_RADIUS + param.sat_altitude
        dist_imt_centre_earth = param.EARTH_RADIUS + param.imt_altitude

        # calculate Cartesian coordinates of satellite, with origin at centre of the Earth
        sat_lat_rad = param.sat_lat_deg * np.pi / 180.
        imt_long_diff_rad = param.imt_long_diff_deg * np.pi / 180.
        x1 = dist_sat_centre_earth * np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth * np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)
        z1 = dist_sat_centre_earth * np.sin(sat_lat_rad)

        # calculate coordinates with origin at IMT system
        imt_lat_rad = param.imt_lat_deg * np.pi / 180.
        satellite_stations.x = [x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)]
        satellite_stations.y = [y1]
        satellite_stations.height = [(z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad)
                                     - dist_imt_centre_earth)]

        satellite_stations.height = [param.sat_altitude]
        satellite_stations.active = True
        satellite_stations.antenna = np.array([AntennaOmni(param.sat_rx_antenna_gain)])
        satellite_stations.bandwidth = param.sat_bandwidth
        satellite_stations.noise_temperature = param.sat_noise_temperature
        satellite_stations.rx_interference = -500

        return satellite_stations
