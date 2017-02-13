# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy
import math

from simulation import Simulation
from parameters import Parameters
from station_manager import StationManager
#from SatelliteReceiver import SatelliteReceiver

class SimulationDownlink(Simulation):
    
    def __init__(self):
        self.imt_base_stations = BaseStations()
#        self.satellite_receiver = SatelliteReceiver()
    
    def initialize(self, *args, **kwargs):
        #self.configure_parameters()
        if Parameters.static_bs:
            self.create_static_bs()
    
    def snapshot(self, *args, **kwargs):
        if not Parameters.static_bs:
            self.create_random_bs()
        self.create_other_station()
        if Parameters.imt_interfered_with:
            self.create_ue()
            self.calculate_path_loss()
            self.connect_ue_to_bs()
            self.apply_power_control()
            self.scheduler()
            self.calculate_powers()
            self.add_external_interference()
            self.recalculate_sinr()
            self.calculate_imt_degradation()
        else:
            #include beamforming
            self.apply_power_control()
            #self.scheduler()
            self.select_active_bs()
            self.calculate_other_interference()
            self.calculate_other_degradation()
        self.collect_results()

    def finalize(self, *args, **kwargs):
        pass
    
    def create_static_bs(self):
        D = Parameters.intersite_distance
        h = D*math.sqrt(3)/2
        
        self.imt_base_stations.x = numpy.array( [ 0, D, D/2, -D/2, -D, -D/2, 
                D/2, 2*D, 3*D/2, D, 0, -D, -3*D/2, -2*D, -3*D/2, -D, 0, D, 3*D/2 ] )
        self.imt_base_stations.y = numpy.array( [ 0, 0, h, h, 0, -h, -h,
                0, h, 2*h, 2*h, 2*h, h, 0, -h, -2*h, -2*h, -2*h, -h ] )
        
    def create_ue(self):
        pass
        
    def select_active_bs(self):
        self.imt_base_stations.active = numpy.random.uniform(size=19) < Parameters.bs_load_probability
    
    def apply_power_control( self ):
        self.imt_base_stations.transmit_power = Parameters.bs_transmit_power
        
    def create_other_station( self ):
        x_min = numpy.min(self.imt_base_stations.x)
        x_max = numpy.max(self.imt_base_stations.x)
        y_min = numpy.min(self.imt_base_stations.y)
        y_max = numpy.max(self.imt_base_stations.y)
                
        self.satellite_receiver.x = ( x_max - x_min ) * numpy.random.random() + x_min
        self.satellite_receiver.y = ( y_max - y_min ) * numpy.random.random() + y_min
        
        