# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
"""

import numpy as np

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
        num_bs = param.num_clusters*param.num_base_stations
        if(param_ant.bs_tx_antenna_type == "BEAMFORMING" or \
           param_ant.bs_rx_antenna_type == "BEAMFORMING"):
            num_bs = 3*num_bs
        imt_base_stations = StationManager(num_bs)
        
        # now we set the coordinates
        if(param_ant.bs_rx_antenna_type == "BEAMFORMING"):
            imt_base_stations.x = np.repeat(topology.x,3)
            imt_base_stations.y = np.repeat(topology.y,3)
        else:
            imt_base_stations.x = topology.x
            imt_base_stations.y = topology.y    
        imt_base_stations.height = param.bs_height*np.ones(num_bs)
        
        imt_base_stations.active = np.ones(num_bs)
        imt_base_stations.tx_power = param.bs_tx_power*np.ones(num_bs)
        imt_base_stations.rx_interference = -500*np.ones(param.ue_k)
        
        if(param_ant.bs_tx_antenna_type == "OMNI"):
            imt_base_stations.tx_antenna = \
                np.array([AntennaOmni(param_ant.bs_tx_omni_antenna_gain) \
                          for i in range(num_bs)])
        elif(param_ant.bs_tx_antenna_type == "BEAMFORMING"):
            imt_base_stations.tx_antenna = np.empty(num_bs,dtype=AntennaBeamformingImt)
            par = param_ant.get_antenna_parameters("BS","TX")
            
            for i in range(num_bs):
                imt_base_stations.tx_antenna[i] = \
                AntennaBeamformingImt(par,param_ant.bs_tx_azimuth[i%3],\
                                      param_ant.bs_tx_elevation)
        
        if(param_ant.bs_rx_antenna_type == "OMNI"):
            imt_base_stations.rx_antenna = \
                np.array([AntennaOmni(param_ant.bs_rx_omni_antenna_gain) \
                          for i in range(num_bs)])
        elif(param_ant.bs_rx_antenna_type == "BEAMFORMING"):
            imt_base_stations.rx_antenna = np.empty(num_bs,dtype=AntennaBeamformingImt)
            par = param_ant.get_antenna_parameters("BS","RX")
            
            for i in range(num_bs):
                imt_base_stations.rx_antenna[i] = \
                AntennaBeamformingImt(par,param_ant.bs_rx_azimuth[i%3],\
                                      param_ant.bs_rx_elevation)
            
        imt_base_stations.bandwidth = param.bandwidth*np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs_noise_figure*np.ones(num_bs)
        imt_base_stations.thermal_noise = -500*np.ones(num_bs)
        imt_base_stations.is_satellite = False
        return imt_base_stations

    @staticmethod
    def generate_imt_ue(param: ParametersImt, 
                        param_ant: ParametersAntennaImt,
                        topology: Topology):
        num_ue = param.num_clusters*param.num_base_stations*param.ue_k*param.ue_k_m
        if(param_ant.bs_tx_antenna_type == "BEAMFORMING" or \
           param_ant.bs_rx_antenna_type == "BEAMFORMING"):
            num_ue = 3*num_ue
        imt_ue = StationManager(num_ue)
        #imt_ue.x = (topology.x_max - topology.x_min)*np.random.random(num_ue) + topology.x_min
        #imt_ue.y = (topology.y_max - topology.y_min)*np.random.random(num_ue) + topology.y_min
        ue_x = list()
        ue_y = list()

        if(param_ant.bs_rx_antenna_type == "BEAMFORMING" or param_ant.bs_rx_antenna_type == "BEAMFORMING"):
            beamforming = True
            cell_r = topology.intersite_distance/3
            num_bs = 3*topology.x.size
            bs_x = np.repeat(topology.x,3)
            bs_y = np.repeat(topology.y,3)
        else:
            beamforming = False
            cell_r = topology.cell_radius
            num_bs = topology.x.size
            bs_x = topology.x
            bs_y = topology.y
            
        for bs in range(num_bs):
            if(beamforming):
                i = bs%3
                ang = np.deg2rad(60+i*120)
                x_center = cell_r*np.cos(ang) + bs_x[bs]
                y_center = cell_r*np.sin(ang) + bs_y[bs]
                # TODO: change this to a parameter that represents the minimum 

                # UE's are generated inside a inscribed circle of a regular hexagon (sector)
                r_max = cell_r*np.sqrt(3)/2
                x_min = x_center - r_max
                x_max = x_center + r_max
                y_min = y_center - r_max
                y_max = y_center + r_max                
                scale = 10
                x = (x_max - x_min)*np.random.random(scale*param.ue_k*param.ue_k_m) + x_min
                y = (y_max - y_min)*np.random.random(scale*param.ue_k*param.ue_k_m) + y_min
            
                r = np.sqrt((x - x_center)**2 + (y - y_center)**2)

                # 2D separation distance between BS and UE
                r_sep_min = param.minimum_separation_distance_bs_ue
                r_sep = np.sqrt((x - bs_x[bs])**2 + (y - bs_y[bs])**2)
                
                idx = np.where((r < r_max) & (r_sep > r_sep_min))[0]
                i = idx[:param.ue_k*param.ue_k_m]
                x = x[i]
                y = y[i]
    
            else:
                x_min = bs_x[bs] - cell_r
                x_max = bs_x[bs] + cell_r
                y_min = bs_y[bs] - cell_r
                y_max = bs_y[bs] + cell_r
                x = (x_max - x_min)*np.random.random(param.ue_k*param.ue_k_m) + x_min
                y = (y_max - y_min)*np.random.random(param.ue_k*param.ue_k_m) + y_min
            ue_x.extend(x)
            ue_y.extend(y)
        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)
        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.height = param.ue_height*np.ones(num_ue)
        imt_ue.tx_power = param.ue_tx_power*np.ones(num_ue)
        imt_ue.rx_interference = -500*np.ones(num_ue)

        if(param_ant.ue_tx_antenna_type == "OMNI"):
            imt_ue.tx_antenna = \
                np.array([AntennaOmni(param_ant.ue_tx_omni_antenna_gain) \
                          for i in range(num_ue)])
        elif(param_ant.ue_tx_antenna_type == "BEAMFORMING"):
            imt_ue.tx_antenna = np.empty(num_ue,dtype=AntennaBeamformingImt)
            par = param_ant.get_antenna_parameters("UE","TX")
            
            if(param_ant.ue_tx_pointing == "FIXED"):
                azi = param_ant.ue_tx_azimuth
                ele = param_ant.ue_tx_elevation
            elif(param_ant.ue_tx_pointing == "RANDOM"):
                azi = np.random.uniform(-180,180)
                ele = np.random.uniform(0,180)
                
            for i in range(num_ue):
                imt_ue.tx_antenna[i] = \
                AntennaBeamformingImt(par,azi,ele)
                
        if(param_ant.ue_rx_antenna_type == "OMNI"):
            imt_ue.rx_antenna = \
                np.array([AntennaOmni(param_ant.ue_rx_omni_antenna_gain) \
                          for i in range(num_ue)])
        elif(param_ant.ue_rx_antenna_type == "BEAMFORMING"):
            imt_ue.rx_antenna = np.empty(num_ue,dtype=AntennaBeamformingImt)
            par = param_ant.get_antenna_parameters("UE","RX")
            
            if(param_ant.ue_rx_pointing == "FIXED"):
                azi = param_ant.ue_rx_azimuth
                ele = param_ant.ue_rx_elevation
            elif(param_ant.ue_rx_pointing == "RANDOM"):
                azi = np.random.uniform(-180,180)
                ele = np.random.uniform(0,180)
                
            for i in range(num_ue):
                imt_ue.rx_antenna[i] = \
                AntennaBeamformingImt(par,azi,ele)
            
        imt_ue.bandwidth = param.bandwidth*np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure*np.ones(num_ue)
        imt_ue.is_satellite = False
        return imt_ue

    @staticmethod
    def generate_fss_stations(param: ParametersFss):
        satellite_stations = StationManager(1)

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
        satellite_stations.rx_antenna = \
            np.array([AntennaOmni(param.sat_rx_antenna_gain)])
        satellite_stations.bandwidth = param.sat_bandwidth
        satellite_stations.noise_temperature = param.sat_noise_temperature
        satellite_stations.rx_interference = -500
        satellite_stations.is_satellite = True

        return satellite_stations
