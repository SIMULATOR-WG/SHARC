# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:03:51 2016

@author: edgar
"""

from sharc.support.observable import Observable
#from support.observer import Observer
from sharc.support.enumerations import State
from sharc.simulation_downlink import SimulationDownlink
from sharc.parameters.parameters_general import ParametersGeneral
from sharc.parameters.parameters_imt import ParametersImt

class Model(Observable):
    """
    Implements the Observable interface. It has a reference to the simulation
    object and controls the simulation flow (init/step/finilize).
    """
    
    def __init__(self):
        super(Model, self).__init__()
        self.simulation = SimulationDownlink(ParametersImt())

    def add_observer(self, observer):
        Observable.add_observer(self, observer)
        self.simulation.add_observer(observer)
        
    def initialize(self):
        """
        Initializes the simulation and performs all pre-simulation tasks
        """
        self.notify_observers(source=__name__,
                              message="Simulation is running...",
                              state=State.RUNNING )
        self.current_snapshot = 1
        self.simulation.initialize()
        
    def step(self):
        """
        Performs one simulation step and collects the results
        """
        if not self.current_snapshot % 10:
            self.notify_observers(source=__name__,
                                  message="Snapshot #" + str(self.current_snapshot))
        #time.sleep(1)
        self.simulation.snapshot()
        self.current_snapshot += 1
            
    def is_finished(self) -> bool:
        """
        Checks is simulation is finished by checking if maximum number of 
        snashots is reached.
        
        Returns
        -------
            True if simulation is finished; False otherwise.
        """
        if self.current_snapshot <= ParametersGeneral.num_snapshots:
            return False
        else:
            return True
            
    def finalize(self):
        """
        Finalizes the simulation and performs all post-simulation tasks
        """
        self.simulation.finalize()
        self.notify_observers(source=__name__, 
                              message="FINISHED!", state=State.FINISHED)
        
    def set_elapsed_time(self, elapsed_time: str):
        """
        Sends the elapsed simulation time to all observers. Simulation time is
        calculated in SimulationThread
        
        Parameters
        ----------
            elapsed_time: Elapsed time.
        """
        self.notify_observers(source=__name__, 
                              message="Elapsed time: " + elapsed_time, 
                              state=State.FINISHED)
