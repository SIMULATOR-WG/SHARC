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

from sharc.support.enumerations import StationType
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

        self.coupling_loss_imt = np.empty(0)
        self.coupling_loss_imt_system = np.empty(0)

        self.bs_to_ue_phi = np.empty(0)
        self.bs_to_ue_theta = np.empty(0)
        self.bs_to_ue_beam_rbs = np.empty(0)

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
        
        self.coupling_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt_system = np.empty(num_ue)

        self.bs_to_ue_phi = np.empty([num_bs, num_ue])
        self.bs_to_ue_theta = np.empty([num_bs, num_ue])
        self.bs_to_ue_beam_rbs = -1.0*np.ones(num_ue, dtype=int)

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
        #self.plot_scenario()
        
        self.connect_ue_to_bs()
        self.select_ue()
        
        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = self.calculate_coupling_loss(self.bs, 
                                                              self.ue,
                                                              self.propagation_imt)
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
            #self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass
        
        self.collect_results(write_to_file, snapshot_number)

        
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
        TODO: calculate coupling loss between activa stations only
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array
        d_2D = station_a.get_distance_to(station_b)
        d_3D = station_a.get_3d_distance_to(station_b)

        if station_b.station_type is StationType.FSS_SS:
            elevation_angles = station_a.get_elevation_angle(station_b, self.param_system)
            path_loss = propagation.get_loss(distance_3D=d_3D, 
                                             frequency=self.param_imt.frequency,
                                             elevation=elevation_angles, 
                                             sat_params = self.param_system,
                                             earth_to_space = True)
        else:
            path_loss = propagation.get_loss(distance_3D=d_3D, 
                                             distance_2D=d_2D, 
                                             frequency=self.param_imt.frequency*np.ones(d_2D.shape),
                                             bs_height=station_a.height,
                                             ue_height=station_b.height,
                                             shadowing=True,
                                             line_of_sight_prob=self.param_imt.line_of_sight_prob)
        # define antenna gains
        gain_a = self.calculate_gains(station_a, station_b)
        gain_b = np.transpose(self.calculate_gains(station_b, station_a))
        
        # calculate coupling loss
        coupling_loss = np.squeeze(path_loss - gain_a - gain_b)
        
        return coupling_loss

        
    def connect_ue_to_bs(self):
        """
        Link the UE's to the serving BS. It is assumed that each group of K*M
        user equipments are distributed and pointed to a certain base station
        according to the decisions taken at TG 5/1 meeting
        """
        num_ue_per_bs = self.param_imt.ue_k*self.param_imt.ue_k_m
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
            self.link[bs] = ue_list


    def select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """               
        self.bs_to_ue_phi, self.bs_to_ue_theta = \
            self.bs.get_pointing_vector_to(self.ue)
        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            # select K UE's among the ones that are connected to BS
            random.shuffle(self.link[bs])
            K = self.param_imt.ue_k
            del self.link[bs][K:]
            # Activate the selected UE's and create beams
            if self.bs.active[bs]:
                self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)
                for ue in self.link[bs]:
                    # add beam to BS antennas
                    self.bs.antenna[bs].add_beam(self.bs_to_ue_phi[bs,ue],
                                             self.bs_to_ue_theta[bs,ue])
                    # add beam to UE antennas
                    self.ue.antenna[ue].add_beam(self.bs_to_ue_phi[bs,ue] - 180,
                                             180 - self.bs_to_ue_theta[bs,ue])
                    # set beam resource block group
                    self.bs_to_ue_beam_rbs[ue] = len(self.bs.antenna[bs].beams_list) - 1

                
    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for
        a given BS
        """
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = self.link[bs]
            self.ue.bandwidth[ue_list] = self.num_rb_per_ue*self.param_imt.rb_bandwidth


    def power_control(self):
        """
        Apply uplink power control algorithm
        """
        if self.param_imt.ue_tx_power_control == "OFF":
            ue_active = np.where(self.ue.active)[0]
            self.ue.tx_power[ue_active] = self.param_imt.ue_tx_power*np.ones(len(ue_active))
        else:
            power_aux =  10*np.log10(self.num_rb_per_ue) + self.param_imt.ue_tx_power_target
            bs_active = np.where(self.bs.active)[0]
            for bs in bs_active:
                ue = self.link[bs]
                power2 = self.coupling_loss_imt[bs,ue]
                self.ue.tx_power[self.link[bs]] = np.minimum(self.param_imt.ue_tx_power,
                                                             self.param_imt.ue_tx_power_alfa*power2 + power_aux)


    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each UE. This is useful only in the
        cases when IMT system is interfered by other system
        """
        # calculate uplink received power for each active BS
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = self.link[bs]
            self.bs.rx_power[bs] = self.ue.tx_power[ue_list] - self.coupling_loss_imt[bs,ue_list]
            # create a list of BSs that serve the interfering UEs
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                ui_list = self.link[bi]
                interference = self.ue.tx_power[ui_list] - self.coupling_loss_imt[bs,ui_list]
                self.bs.rx_interference[bs] = 10*np.log10( \
                    np.power(10, 0.1*self.bs.rx_interference[bs])
                    + np.power(10, 0.1*interference))

            # calculate N
            self.bs.thermal_noise[bs] = \
                10*np.log10(self.param_imt.BOLTZMANN_CONSTANT*self.param_imt.noise_temperature) + \
                10*np.log10(self.num_rb_per_ue*self.param_imt.rb_bandwidth * 1e6) + \
                self.bs.noise_figure[bs]
    
            # calculate I+N
            self.bs.total_interference[bs] = \
                10*np.log10(np.power(10, 0.1*self.bs.rx_interference[bs]) + \
                            np.power(10, 0.1*self.bs.thermal_noise[bs]))
                
            # calculate SNR and SINR
            self.bs.sinr[bs] = self.bs.rx_power[bs] - self.bs.total_interference[bs]
            self.bs.snr[bs] = self.bs.rx_power[bs] - self.bs.thermal_noise[bs]


    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """

        self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system, 
                                                                     self.ue,
                                                                     self.propagation_system)

        ue_bandwidth = self.num_rb_per_ue * self.param_imt.rb_bandwidth

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        ue_active = np.where(self.ue.active)[0]
        interference_ue = self.ue.tx_power[ue_active] - self.coupling_loss_imt_system[ue_active] \
                            + 10*math.log10(ue_bandwidth/self.param_system.sat_bandwidth)

        # calculate the aggregate interference on system
        self.system.rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference_ue)))

        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.sat_noise_temperature) + \
                          10*math.log10(self.param_system.sat_bandwidth * 1e6)

        # calculate I+N
        self.system.total_interference = \
            10*np.log10(np.power(10, 0.1*self.system.rx_interference) + \
                        np.power(10, 0.1*self.system.thermal_noise))

        # calculate INR at the system
        self.system.inr = self.system.rx_interference - self.system.thermal_noise


    def collect_results(self, write_to_file: bool, snapshot_number: int):
        self.results.imt_ul_coupling_loss.extend( \
            np.reshape(self.coupling_loss_imt, self.ue.num_stations*self.bs.num_stations).tolist())
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
        phi, theta = station_a.get_pointing_vector_to(station_b)
        
        if(station_a.station_type == StationType.IMT_BS):
            beams_idx = self.bs_to_ue_beam_rbs
        elif(station_a.station_type == StationType.IMT_UE):
            beams_idx = np.zeros(self.bs.num_stations,dtype=int)
        elif(station_a.station_type == StationType.FSS_SS):
            beams_idx = np.zeros(self.ue.num_stations,dtype=int)
        
        gains = np.zeros(phi.shape)
        station_a_active = np.where(station_a.active)[0]
        station_b_active = np.where(station_b.active)[0]
        for k in station_a_active:
            gains[k,station_b_active] = station_a.antenna[k].calculate_gain(phi_vec=phi[k,station_b_active],
                                                                            theta_vec=theta[k,station_b_active],
                                                                            beams_l=beams_idx[station_b_active])
                
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
        
        
    def plot_scenario(self):
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        
        # Plot network topology
        self.topology.plot(ax)
        
        # Plot user equipments
        ax.scatter(self.ue.x, self.ue.y, color='r', edgecolor="w", linewidth=0.5, label="UE")
        
        # Plot UE's azimuth
        d = 0.1 * self.topology.cell_radius
        for i in range(len(self.ue.x)):
            plt.plot([self.ue.x[i], self.ue.x[i] + d*math.cos(math.radians(self.ue.azimuth[i]))], 
                     [self.ue.y[i], self.ue.y[i] + d*math.sin(math.radians(self.ue.azimuth[i]))], 
                     'r-')        
        
        plt.axis('image') 
        plt.title("Simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()    
        plt.show()        
        
        sys.exit(0)
        