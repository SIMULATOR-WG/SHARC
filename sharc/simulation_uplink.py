# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
"""

import numpy as np
import random
import math
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_factory import StationFactory
from sharc.station_manager import StationManager
from sharc.topology.topology_factory import TopologyFactory
from sharc.propagation.propagation_factory import PropagationFactory
from sharc.propagation.propagation import Propagation
from sharc.results import Results

class SimulationUplink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param_imt: ParametersImt, param_system: ParametersFss, param_ant: ParametersAntennaImt):
        super().__init__()
        self.param_imt = param_imt
        self.param_system = param_system
        self.param_imt_antenna = param_ant

        self.topology = TopologyFactory.createTopology(self.param_imt)

        self.propagation_imt = PropagationFactory.createPropagation(self.param_imt.channel_model)
        self.propagation_system = PropagationFactory.createPropagation(self.param_system.channel_model)

        self.path_loss_imt = np.empty(0)
        self.path_loss_imt_system = np.empty(0)
        self.coupling_loss_imt = np.empty(0)
        self.coupling_loss_imt_system = np.empty(0)

        self.phi = np.empty(0)
        self.theta = np.empty(0)

        self.ue = np.empty(0)
        self.bs = np.empty(0)
        self.system = np.empty(0)
        
        self.link = dict()

        self.num_rb_per_bs = 0
        self.num_rb_per_ue = 0

        self.results = Results()

        
    def initialize(self, *args, **kwargs):
        self.topology.calculate_coordinates()
        num_bs = self.topology.num_base_stations
        num_ue = num_bs*self.param_imt.ue_k*self.param_imt.ue_k_m
        
        self.path_loss_imt = np.empty([num_bs, num_ue])
        self.path_loss_imt_system = np.empty(num_ue)
        self.coupling_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt_system = np.empty(num_ue)

        self.phi = np.empty([num_bs, num_ue])
        self.theta = np.empty([num_bs, num_ue])

        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.system = np.empty(1)

        # this attribute indicates the list of UE's that are connected to each
        # base station. The position the the list indicates the resource block
        # group that is allocated to the given UE
        self.link = dict([(bs,list()) for bs in range(num_bs)])

        # calculates the number of RB per BS
        self.num_rb_per_bs = math.trunc((1-self.param_imt.guard_band_ratio)* \
                            self.param_imt.bandwidth /self.param_imt.rb_bandwidth)
        # calculates the number of RB per UE on a given BS
        self.num_rb_per_ue = math.trunc(self.num_rb_per_bs/self.param_imt.ue_k)        
        
        
    def snapshot(self, write_to_file, snapshot_number, *args, **kwargs):
        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates()
        
        # Create the base stations (remember that it takes into account the
        # network load factor)
        self.bs = StationFactory.generate_imt_base_stations(self.param_imt,
                                                            self.param_imt_antenna,
                                                            self.topology)      
        
        # Create the other system (FSS, HAPS, etc...)
        # Currently it supports only FSS space station
        self.system = StationFactory.generate_fss_stations(self.param_system)

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.param_imt,
                                                 self.param_imt_antenna,
                                                 self.topology)
        #self.plot_macrocell_scenario()
        #self.plot_hotspot_scenario()
        #sys.exit(0)
        
        # reset the index of beams
        #self.beams_idx = -1*np.ones(self.ue.num_stations, dtype=int)
        
        self.connect_ue_to_bs()
        self.select_ue()
        
        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = np.transpose( \
                             self.calculate_coupling_loss(self.ue, self.bs,
                                                          self.propagation_imt))
        self.scheduler()
        self.power_control()
        
        if self.param_imt.interfered_with:
            # Execute this piece of code if the other system generates 
            # interference into IMT
            #self.calculate_sinr()
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
            pass
        else:
            # Execute this piece of code if IMT generates interference into
            # the other system
            self.calculate_sinr()
            self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass
        
        self.collect_results(write_to_file, snapshot_number)
        self.reset_antennas()

        
    def finalize(self, snapshot_number, *args, **kwargs):
        self.results.write_files(snapshot_number)
        self.notify_observers(source=__name__, results=self.results)


    def calculate_coupling_loss(self,
                                station_a: StationManager,
                                station_b: StationManager,
                                propagation: Propagation) -> np.array:
        """
        Calculates the path coupling loss from each station_a to all station_b.
        Result is returned as a numpy array with dimensions num_a x num_b
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array
        d_2D = station_a.get_distance_to(station_b)
        d_3D = station_a.get_3d_distance_to(station_b)

        if station_b.is_satellite:
            elevation_angles = station_a.get_elevation_angle(station_b, self.param_system)
            self.path_loss = propagation.get_loss(distance=d_3D, frequency=self.param_imt.frequency,
                                             elevation=elevation_angles, sat_params = self.param_system,
                                             earth_to_space = True)
        else:
            self.path_loss = propagation.get_loss(distance=np.transpose(d_3D), 
                                                  distance_2D=np.transpose(d_2D), 
                                                  frequency=self.param_imt.frequency*np.ones(np.transpose(d_3D).shape),
                                                  bs_height=station_b.height,
                                                  ue_height=station_a.height,
                                                  shadowing=False)
        # define antenna gains
        gain_a = self.calculate_gains(station_a,station_b,"TX")
        gain_b = self.calculate_gains(station_b,station_a,"RX")
        if(station_b.num_stations > 1):
            gain_b_t = np.transpose(gain_b)
        else:
            gain_b_t = gain_b[:,np.newaxis]
        # calculate coupling loss
        return self.path_loss - gain_a - gain_b_t

        
    def connect_ue_to_bs(self):
        """
        Link the UE's to the serving BS. It is assumed that each group of K*M
        user equipments are distributed and pointed to a certain base station
        according to the decisions taken at TG 5/1 meeting
        """
        num_ue_per_bs = self.param_imt.ue_k*self.param_imt.ue_m
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
            self.link[bs] = ue_list


    def select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            # select K UE's among the ones that are connected to BS
            random.shuffle(self.link[bs])
            K = self.param_imt.ue_k
            del self.link[bs][K:]
            # Activate the selected UE's
            if self.bs.active[bs]:
                self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)
            for ue in self.link[bs]:
                # add beam to antennas
                self.ue.antenna[ue].add_beam(self.phi[bs,ue] - 180,
                                             180 - self.theta[bs,ue])

                
    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for
        a given BS
        """
        for bs in range(self.bs.num_stations):
            ue_list = self.link[bs]
            self.ue.bandwidth[ue_list] = self.num_rb_per_ue*self.param_imt.rb_bandwidth


    def power_control(self):
        """
        Apply uplink power control algorithm
        """
        if self.param_imt.ue_tx_power_control == "OFF":
            self.param_imt.ue.tx_power = self.param_imt.ue_tx_power*np.ones(self.ue.num_stations)
        else:
            power_aux =  10*np.log10(self.num_rb_per_ue) + self.param_imt.ue_tx_power_target
            for bs in range(self.bs.num_stations):
                power2 = self.coupling_loss[self.link[bs], bs]
                self.ue.tx_power[self.link[bs]] = np.minimum(self.param_imt.ue_tx_power,
                                                             self.param_imt.ue_tx_power_alfa*power2 + power_aux)


    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each UE. This is useful only in the
        cases when IMT system is interfered by other system
        """
        #bs_all = [b for b in range(self.bs.num_stations)]
        bs_active = np.where(self.bs.active)[0]

        self.bs.rx_power = dict([(bs, -500 * np.ones(len(self.link[bs]))) for bs in bs_active])
        self.bs.rx_interference = dict([(bs, -500 * np.ones(len(self.link[bs]))) for bs in bs_active])
        self.bs.total_interference = dict([(bs, -500 * np.ones(len(self.link[bs]))) for bs in bs_active])
        self.bs.snr = dict([(bs, -500 * np.ones(len(self.link[bs]))) for bs in bs_active])
        self.bs.sinr = dict([(bs, -500 * np.ones(len(self.link[bs]))) for bs in bs_active])

        # calculate uplink received power for each active BS
        for bs in bs_active:
            ue_list = self.link[bs]
            self.bs.rx_power[bs] = self.ue.tx_power[ue_list] - self.coupling_loss[bs,ue_list]
            # create a list of BSs that serve the interfering UEs
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                ui_list = self.link[bi]
                interference = self.ue.tx_power[ui_list] - self.coupling_loss[bs,ui_list]
                self.bs.rx_interference[bs] = 10*np.log10( \
                    np.power(10, 0.1*self.bs.rx_interference[bs])
                    + np.power(10, 0.1*interference))

            self.bs.thermal_noise[bs] = \
                10*np.log10(self.param_imt.BOLTZMANN_CONSTANT*self.param_imt.noise_temperature) + \
                10*np.log10(self.num_rb_per_ue*self.param_imt.rb_bandwidth * 1e6) + \
                self.bs.noise_figure[bs]
    
            self.bs.total_interference[bs] = \
                10*np.log10(np.power(10, 0.1*self.bs.rx_interference[bs]) + \
                            np.power(10, 0.1*self.bs.thermal_noise[bs]))
                
            self.bs.sinr[bs] = self.bs.rx_power[bs] - self.bs.total_interference[bs]
            self.bs.snr[bs] = self.bs.rx_power[bs] - self.bs.thermal_noise[bs]


    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """

        self.coupling_loss_ue_sat = np.array(np.transpose(
                                self.calculate_coupling_loss(self.ue, self.system,
                                            self.propagation_system)).tolist()[0])
        self.beams_idx = -1*np.ones(self.ue.num_stations,dtype=int)

        ue_bandwidth = self.num_rb_per_ue * self.param_imt.rb_bandwidth

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        ue_active = np.where(self.ue.active)[0]
        interference_ue = self.ue.tx_power[ue_active] - self.coupling_loss_ue_sat[ue_active] \
                            + 10*math.log10(ue_bandwidth/self.param_system.sat_bandwidth)

        self.system.rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference_ue)))

        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.sat_noise_temperature) + \
                          10*math.log10(self.param_system.sat_bandwidth * 1e6)

        self.system.total_interference = \
            10*np.log10(np.power(10, 0.1*self.system.rx_interference) + \
                        np.power(10, 0.1*self.system.thermal_noise))

        self.system.inr = self.system.rx_interference - self.system.thermal_noise


    def collect_results(self, write_to_file: bool, snapshot_number: int):
        self.results.imt_ul_coupling_loss.extend( \
            np.reshape(self.coupling_loss, self.ue.num_stations*self.bs.num_stations).tolist())
        self.results.system_inr.extend([self.system.inr])
        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = self.link[bs]
            tput = self.calculate_imt_ul_tput(self.bs.sinr[bs])
            self.results.imt_ul_tput.extend(tput.tolist())
            self.results.imt_ul_tx_power.extend(self.ue.tx_power[ue_list].tolist())
            imt_ul_tx_power_density = 10*np.log10(np.power(10, 0.1*self.ue.tx_power[ue_list])/(self.num_rb_per_ue*self.param_imt.rb_bandwidth*1e6))
            self.results.imt_ul_tx_power_density.extend(imt_ul_tx_power_density.tolist())
            self.results.imt_ul_sinr.extend(self.bs.sinr[bs].tolist())
            self.results.imt_ul_snr.extend(self.bs.snr[bs].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)

            
    def calculate_gains(self,
                        station_a: StationManager,
                        station_b: StationManager) -> np.array:
        """
        Calculates the gains of antennas in station_a in the direction of
        station_b       
        """
        if(station_a.num_stations > 1):
            point_vec_x = station_b.x- station_a.x[:,np.newaxis]
            point_vec_y = station_b.y - station_a.y[:,np.newaxis]
            point_vec_z = station_b.height - station_a.height[:,np.newaxis]
        else:
            point_vec_x = station_b.x- station_a.x
            point_vec_y = station_b.y - station_a.y
            point_vec_z = station_b.height - station_a.height
            
        dist = station_a.get_3d_distance_to(station_b)
        
        self.phi = np.rad2deg(np.arctan2(point_vec_y,point_vec_x))
        self.theta = np.rad2deg(np.arccos(point_vec_z/dist))
        
        gains = np.zeros_like(self.phi)
        if(len(np.shape(gains)) != 1):
            for k in range(station_a.num_stations):
                gains[k,:] = station_a.antenna[k].calculate_gain(self.phi[k,:],\
                     self.theta[k,:],self.beams_idx)
        else:
            gains = station_a.tx_antenna[0].calculate_gain(self.phi,self.theta)
                
        return gains
    
        
    def calculate_imt_ul_tput(self, sinr: np.array) -> np.array:
        tput_min = 0
        tput_max = self.param_imt.ul_attenuation_factor*math.log2(1+math.pow(10, 0.1*self.param_imt.ul_sinr_max))
        
        tput = self.param_imt.ul_attenuation_factor*np.log2(1+np.power(10, 0.1*sinr))
        
        id_min = np.where(sinr<self.param_imt.ul_sinr_min)[0]
        id_max = np.where(sinr>self.param_imt.ul_sinr_max)[0]

        if len(id_min) > 0:
            tput[id_min] = tput_min
        if len(id_max) > 0:
            tput[id_max] = tput_max

        return tput
        
        
    def plot_macrocell_scenario(self):
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        
        #plot hexagons
        r = self.topology.intersite_distance/3
        for x, y, az in zip(self.topology.x, self.topology.y, self.topology.azimuth):
            se = list([[x,y]])
            angle = int(az - 60)
            for a in range(6):
                se.extend([[se[-1][0] + r*math.cos(math.radians(angle)), se[-1][1] + r*math.sin(math.radians(angle))]])
                angle += 60
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)
        
        # macro cell base stations
        plt.scatter(self.topology.x, self.topology.y, color='k', edgecolor="k", linewidth=4, label="BS")
        
        # UE's
        plt.scatter(self.ue.x, self.ue.y, color='r', edgecolor="w", linewidth=0.5, label="UE")
        
#        # UE azimuth
#        d = 0.2 * self.topology.cell_radius
#        for i in range(len(self.ue.x)):
#            plt.plot([self.ue.x[i], self.ue.x[i] + d*math.cos(math.radians(self.ue.azimuth[i]))], 
#                     [self.ue.y[i], self.ue.y[i] + d*math.sin(math.radians(self.ue.azimuth[i]))], 
#                     'r-')
        
        # plot macro cell coverage area
    #    r = (topology.macrocell.intersite_distance/3)*math.sqrt(3)/2 - topology.param.max_dist_hotspot_ue/2
    #    for x, y, az in zip(topology.macrocell.x, topology.macrocell.y, topology.macrocell.azimuth):
    #        # find the center coordinates of the sector (hexagon)
    #        mx = x + topology.macrocell.intersite_distance/3*math.cos(math.radians(az))
    #        my = y + topology.macrocell.intersite_distance/3*math.sin(math.radians(az))
    #        circ = plt.Circle((mx, my), radius=r, color='b', fill=False, linewidth=0.5)
    #        ax.add_patch(circ)    
        
        plt.axis('image') 
        plt.title("Simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        #plt.xlim((-3000, 3000))
        #plt.ylim((-3000, 3000))                
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()    
        plt.show()
            
        
    def plot_hotspot_scenario(self):
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        
        #plot hexagons
        r = self.topology.macrocell.intersite_distance/3
        for x, y, az in zip(self.topology.macrocell.x, self.topology.macrocell.y, self.topology.macrocell.azimuth):
            se = list([[x,y]])
            angle = int(az - 60)
            for a in range(6):
                se.extend([[se[-1][0] + r*math.cos(math.radians(angle)), se[-1][1] + r*math.sin(math.radians(angle))]])
                angle += 60
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)
        
        # macro cell base stations
        plt.scatter(self.topology.macrocell.x, self.topology.macrocell.y, color='k', edgecolor="k", linewidth=4, label="BS")
        
        # plot hotspots
        plt.scatter(self.topology.x, self.topology.y, color='g', edgecolor="w", linewidth=0.5, label="Hotspot")        
        
        # UE's
        plt.scatter(self.ue.x, self.ue.y, color='r', edgecolor="w", linewidth=0.5, label="UE")
        
#        # UE azimuth
#        d = 0.2 * self.topology.cell_radius
#        for i in range(len(self.ue.x)):
#            plt.plot([self.ue.x[i], self.ue.x[i] + d*math.cos(math.radians(self.ue.azimuth[i]))], 
#                     [self.ue.y[i], self.ue.y[i] + d*math.sin(math.radians(self.ue.azimuth[i]))], 
#                     'r-')
        
        # plot hotspots coverage area
        for x, y, a in zip(self.topology.x, self.topology.y, self.topology.azimuth):
            pa = patches.Wedge( (x, y), self.topology.cell_radius, a-60, a+60, fill=False, 
                               edgecolor="green", linestyle='solid' )
            ax.add_patch(pa)        
        
        plt.axis('image') 
        plt.title("Hotspots simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        #plt.xlim((-3000, 3000))
        #plt.ylim((-3000, 3000))                
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()    
        plt.show()        