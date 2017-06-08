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

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_factory import StationFactory
from sharc.station_manager import StationManager
from sharc.topology.topology_factory import TopologyFactory
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_close_in import PropagationCloseIn
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.propagation.propagation_sat_simple import PropagationSatSimple
from sharc.propagation.propagation_uma import PropagationUMa


from sharc.results import Results

class SimulationUplink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param: ParametersImt, param_system: ParametersFss, param_ant: ParametersAntennaImt):
        super(SimulationUplink, self).__init__()
        self.param = param
        self.param_imt_antenna = param_ant
        self.param_system = param_system

        self.topology = TopologyFactory.createTopology(self.param)

        if self.param.channel_model == "FSPL":
            self.propagation_imt = PropagationFreeSpace()
        elif self.param.channel_model == "UMa":
            self.propagation_imt = PropagationUMa()
        elif self.param.channel_model == "CI":
            self.propagation_imt = PropagationCloseIn(self.param.topology,
                                                      self.param.line_of_sight_prob)
        else:
            sys.stderr.write("error: invalid parameter channel_model" + self.param.channel_model
                             + "in IMT propagation\n")
            sys.exit(1)


        if self.param_system.channel_model == "FSPL":
            self.propagation_system = PropagationFreeSpace()
        elif self.param_system.channel_model == "SatelliteSimple":
            self.propagation_system = PropagationSatSimple(self.param_system.line_of_sight_prob)
        elif self.param_system.channel_model == "P619":
            self.propagation_system = PropagationP619()
        else:
            sys.stderr.write("error: invalid parameter channel_model" + self.param_system.channel_model
                             + "in satellite propagation\n")
            sys.exit(1)

        num_ue = self.param.num_clusters*self.param.num_base_stations \
                                 *self.param.ue_k*self.param.ue_k_m
        num_bs = self.param.num_clusters*self.param.num_base_stations
        
        # Emulate 3 cells per site by multiplying the number of BSs by 3
        if(self.param_imt_antenna.bs_rx_antenna_type == "BEAMFORMING"):
            num_bs = 3*num_bs
        if(self.param_imt_antenna.ue_tx_antenna_type == "BEAMFORMING"):
            num_ue = 3*num_ue

        self.interference_ue = np.empty(num_ue)
        self.coupling_loss = np.empty([num_bs, num_ue])
        self.coupling_loss_ue_sat = np.empty(num_ue)
        self.coupling_loss_bs_sat = np.empty(num_bs)

        self.phi = np.empty([num_bs, num_ue])
        self.theta = np.empty([num_bs, num_ue])
        self.path_loss = np.empty(num_ue)

        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.bs_load_prob = param.bs_load_probability
        self.system = np.empty(1)

        self.ue_tx_power = np.empty([num_ue, num_bs])
        self.ue_tx_power_control = param.ue_tx_power_control
        self.ue_tx_power_target = param.ue_tx_power_target
        self.ue_tx_power_alfa = param. ue_tx_power_alfa
        # this attribute indicates the list of UE's that are connected to each
        # base station. The position the the list indicates the resource block
        # group that is allocated to the given UE
        self.link = dict([(bs,list()) for bs in range(num_bs)])

        # calculates the number of RB per BS
        self.num_rb_per_bs = math.trunc((1-self.param.guard_band_ratio)* \
                            self.param.bandwidth /self.param.rb_bandwidth)
        # calculates the number of RB per UE on a given BS
        self.num_rb_per_ue = math.trunc(self.num_rb_per_bs/self.param.ue_k)

        self.results = Results()

    def initialize(self, *args, **kwargs):
        self.bs = StationFactory.generate_imt_base_stations(self.param,
                                                            self.param_imt_antenna,
                                                            self.topology)

    def snapshot(self, write_to_file, snapshot_number, *args, **kwargs):
        self.create_system()
        self.update_bs()
        self.create_ue()
        self.coupling_loss = np.transpose( \
                             self.calculate_coupling_loss(self.ue, self.bs,
                                                          self.propagation_imt))
        self.connect_ue_to_bs()
        self.select_ue()
        self.scheduler()
        self.power_control()

        if not self.param.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass

        if self.param.interfered_with:
            pass
            #self.calculate_sinr()
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
        else:
            self.calculate_sinr()
            self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass
        
        self.collect_results(write_to_file, snapshot_number)
        self.reset_antennas()

    def finalize(self, snapshot_number, *args, **kwargs):
        self.results.write_files(snapshot_number)
        self.notify_observers(source=__name__, results=self.results)

    def create_ue(self):
        self.ue = StationFactory.generate_imt_ue(self.param, \
                                                 self.param_imt_antenna,\
                                                 self.topology)
        #self.plot_scenario(self.topology.intersite_distance/3, self.bs.x, self.bs.y, self.ue.x, self.ue.y)
        #sys.exit(1)

    def update_bs(self):
        self.bs.active = np.random.rand(self.bs.num_stations) < self.bs_load_prob


    def create_system(self):
        self.system = StationFactory.generate_fss_stations(self.param_system)

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
            self.path_loss = propagation.get_loss(distance=d_3D, frequency=self.param.frequency,
                                             elevation=elevation_angles, sat_params = self.param_system,
                                             earth_to_space = True)
        else:
            self.path_loss = propagation.get_loss(distance=np.transpose(d_3D), 
                                                  distance_2D=np.transpose(d_2D), 
                                                  frequency=self.param.frequency*np.ones(np.transpose(d_3D).shape),
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
#        self.coupling_loss = np.maximum(path_loss - tx_gain - rx_gain,
#                                          ParametersImt.mcl)

    def connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within
        the minimum coupling loss plus the HO margin.
        """
        for ue in range(self.ue.num_stations):
            # for each receiver, select the base stations
            minimum_coupling_loss = np.amin(self.coupling_loss[:,ue])
            bs_all = np.where(self.coupling_loss[:,ue] < minimum_coupling_loss + self.param.ho_margin)[0]
            bs = np.random.choice(bs_all)
            self.link[bs].append(ue)
            # add beam to antennas
            if(self.param_imt_antenna.bs_rx_antenna_type == "BEAMFORMING"):
                self.bs.rx_antenna[bs].add_beam(self.phi[bs,ue],self.theta[bs,ue])
            if(self.param_imt_antenna.ue_tx_antenna_type == "BEAMFORMING"):
                self.ue.tx_antenna[ue].add_beam(self.phi[bs,ue] - 180,\
                                  180 - self.theta[bs,ue])

    def select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        for bs in range(self.bs.num_stations):
            # select K UE's among the ones that are connected to BS
            random.shuffle(self.link[bs])
            K = self.param.ue_k
            #K = np.random.choice(range(1, len(self.link[bs])))
            del self.link[bs][K:]
            # Activate the selected UE's
            if self.bs.active[bs]:
                self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)

    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for
        a given BS
        """
        for bs in range(self.bs.num_stations):
            ue_list = self.link[bs]
            self.ue.bandwidth[ue_list] = self.num_rb_per_ue*self.param.rb_bandwidth

    def power_control(self):
        """
        Apply uplink power control algorithm
        """

        if self.ue_tx_power_control == "OFF":
            self.ue.tx_power = self.param.ue_tx_power*np.ones(self.ue.num_stations)
        else:
            power_aux =  10*np.log10(self.num_rb_per_ue) + self.ue_tx_power_target
            for bs in range(self.bs.num_stations):
                power2 = self.path_loss[self.link[bs], bs]
                self.ue.tx_power[self.link[bs]] = np.minimum(self.param.ue_tx_power,\
                self.ue_tx_power_alfa*power2+power_aux)

                power2 = self.path_loss[self.link[bs], bs]
                self.ue.tx_power[self.link[bs]] = np.minimum(self.param.ue_tx_power,\
                self.ue_tx_power_alfa*power2+power_aux)

    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each UE. This is useful only in the
        cases when IMT system is interfered by other system
        """
        #bs_all = [b for b in range(self.bs.num_stations)]
        bs_active = np.where( self.bs.active )[0]

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
                10*np.log10(self.param.BOLTZMANN_CONSTANT*self.param.noise_temperature) + \
                10*np.log10(self.num_rb_per_ue*self.param.rb_bandwidth * 1e6) + \
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

        ue_bandwidth = self.num_rb_per_ue * self.param.rb_bandwidth

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
            imt_ul_tx_power_density = 10*np.log10(np.power(10, 0.1*self.ue.tx_power[ue_list])/(self.num_rb_per_ue*self.param.rb_bandwidth*1e6))
            self.results.imt_ul_tx_power_density.extend(imt_ul_tx_power_density.tolist())
            self.results.imt_ul_sinr.extend(self.bs.sinr[bs].tolist())
            self.results.imt_ul_snr.extend(self.bs.snr[bs].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)

    def calculate_gains(self,
                        station_a: StationManager,
                        station_b: StationManager,
                        antenna_txrx: str) -> np.array:
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
        if(antenna_txrx == "TX"):
            if(len(np.shape(gains)) != 1):
                for k in range(station_a.num_stations):
                    gains[k,:] = station_a.tx_antenna[k].calculate_gain(self.phi[k,:],\
                         self.theta[k,:])
            else:
                gains = station_a.tx_antenna[0].calculate_gain(self.phi,self.theta)
        elif(antenna_txrx == "RX"):
            if(len(np.shape(gains)) != 1):
                for k in range(station_a.num_stations):
                    gains[k,:] = station_a.rx_antenna[k].calculate_gain(self.phi[k,:],\
                         self.theta[k,:])
            else:
                gains = station_a.rx_antenna[0].calculate_gain(self.phi,self.theta)
                
        return gains
    
    def reset_antennas(self):
        if(self.param_imt_antenna.bs_rx_antenna_type == "BEAMFORMING"):
            for bs in range(self.bs.num_stations):
                self.bs.rx_antenna[bs].reset_beams()
        if(self.param_imt_antenna.ue_tx_antenna_type == "BEAMFORMING"):
            for ue in range(self.ue.num_stations):
                self.ue.tx_antenna[ue].reset_beams()

    def calculate_imt_ul_tput(self, sinr: np.array) -> np.array:
        tput_min = 0
        tput_max = self.param.ul_attenuation_factor*math.log2(1+math.pow(10, 0.1*self.param.ul_sinr_max))
        
        tput = self.param.ul_attenuation_factor*np.log2(1+np.power(10, 0.1*sinr))
        
        id_min = np.where(sinr<self.param.ul_sinr_min)[0]
        id_max = np.where(sinr>self.param.ul_sinr_max)[0]

        if len(id_min) > 0:
            tput[id_min] = tput_min
        if len(id_max) > 0:
            tput[id_max] = tput_max

        return tput
        
    def plot_scenario(self, cell_radius, bs_x, bs_y, ue_x, ue_y):
        psi = np.radians([60, 120, 240, 300])
    
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        
        r = cell_radius
        for x, y in zip(bs_x, bs_y):
            se = list([[x,y]])
            se.extend([[se[-1][0] + r, se[-1][1]]])
            se.extend([[se[-1][0] + r*math.cos(psi[0]), se[-1][1] + r*math.sin(psi[0])]])
            se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
            se.extend([[se[-1][0] - r, se[-1][1]]])
            se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)
    
            se = list([[x,y]])
            se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
            se.extend([[se[-1][0] - r, se[-1][1]]])
            se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
            se.extend([[se[-1][0] + r*math.cos(psi[3]), se[-1][1] + r*math.sin(psi[3])]])
            se.extend([[se[-1][0] + r, se[-1][1]]])
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)
            
            se = list([[x,y]])
            se.extend([[se[-1][0] + r, se[-1][1]]])
            se.extend([[se[-1][0] + r*math.cos(psi[3]), se[-1][1] + r*math.sin(psi[3])]])
            se.extend([[se[-1][0] + r*math.cos(psi[2]), se[-1][1] + r*math.sin(psi[2])]])
            se.extend([[se[-1][0] - r, se[-1][1]]])
            se.extend([[se[-1][0] + r*math.cos(psi[1]), se[-1][1] + r*math.sin(psi[1])]])
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)       
            
    
        # plot hotspot centers
        #plt.scatter(topology.hotspot_x, topology.hotspot_y, color='k', edgecolor="w", linewidth=0.5)
        
        # plot small cells
        #plt.scatter(topology.x, topology.y, color='r', edgecolor="w", linewidth=0.5, label="Small cell")
        
        # plot hotspots coverage area
#        for hx, hy in zip(topology.hotspot_x, topology.hotspot_y):
#            circ = plt.Circle((hx, hy), radius=50, color='g', fill=False, linewidth=0.5)
#            ax.add_patch(circ)
        
        # macro cell base stations
        plt.scatter(bs_x, bs_y, color='k', edgecolor="k", linewidth=4, label="BS")

        # UE
        plt.scatter(ue_x, ue_y, color='r', edgecolor="w", linewidth=0.5, label="UE")
        
        # sector centers
        #plt.scatter(-sector_y, sector_x, color='g', edgecolor="g")
        
        # plot macro cell coverage area
        #ax = fig.gca()
    #    for mx, my in zip(topology.topology_macrocell.x, topology.topology_macrocell.y):
    #        circ = plt.Circle((mx, my), radius=666.667*math.sqrt(3)/2-70, color='b', fill=False, linewidth=0.5)
    #        ax.add_patch(circ)  
    
        # plot separation radius
        for mx, my in zip(bs_x, bs_y):
            circ = plt.Circle((mx, my), radius=10, color='g', fill=False, linewidth=0.5)
            ax.add_patch(circ)  
        
    
    
        
        plt.axis('image') 
        plt.title("Macro cell topology")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        #plt.xlim((-3000, 3000))
        #plt.ylim((-3000, 3000))                
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()    
        plt.show()
            