# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:03:51 2016

@author: edgar
"""

import time

from support.observable import Observable
from support.enumerations import State
#from simulation_downlink import SimulationDownlink
from parameters import Parameters

class Model(Observable):
    """
    Implements the Observable interface. It has a reference to the simulation
    object and controls the simulation flow (init/step/finilize).
    
    Attributes
    ----------
        __simulation : Simulation
    """
    
    def __init__(self):
        super(Model, self).__init__()
        #self.simulation = SimulationDownlink()
        
    def initialize(self):
        self.notify_observers(source=__name__,
                              message="Simulation is running...",
                              state=State.RUNNING )
        self.current_snapshot = 1
        #self.simulation.initialize()
        
    def step(self):
        self.notify_observers(source=__name__,
                              message="Snapshot #" + str(self.current_snapshot))
        time.sleep(1)
        #self.simulation.snapshot()
        self.current_snapshot += 1
            
    def is_finished(self):
        if self.current_snapshot <= Parameters.num_snapshots:
            return False
        else:
            return True
            
    def finalize(self):
        #self.simulation.finalize()
        self.notify_observers(source=__name__, 
                              message="FINISHED!", state=State.FINISHED)
        
    def set_elapsed_time(self, elapsed_time):
        self.notify_observers(source=__name__, 
                              message="Elapsed time: " + elapsed_time, 
                              state=State.FINISHED)
