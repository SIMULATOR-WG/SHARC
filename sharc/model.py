# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:03:51 2016

@author: edgar
"""

from sharc.support.observable import Observable
#from support.observer import Observer
from sharc.support.enumerations import State
from sharc.simulation_downlink import SimulationDownlink
from sharc.simulation_uplink import SimulationUplink
from sharc.parameters.parameters_general import ParametersGeneral
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.parameters.parameters_fss_es import ParametersFssEs

class Model(Observable):
    """
    Implements the Observable interface. It has a reference to the simulation
    object and controls the simulation flow (init/step/finilize).
    """
    
    def __init__(self):
        super(Model, self).__init__()
        if ParametersGeneral().system == "FSS_SS":
            param_system = ParametersFss()
        if ParametersGeneral().system == "FSS_ES":
            param_system = ParametersFssEs()
        
        if ParametersGeneral.imt_link == "DOWNLINK":
            self.simulation = SimulationDownlink(ParametersImt(), param_system, ParametersAntennaImt())
        else:
            self.simulation = SimulationUplink(ParametersImt(), param_system, ParametersAntennaImt())

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
        self.current_snapshot = 0
        self.simulation.initialize()
        
    def snapshot(self):
        """
        Performs one simulation step and collects the results
        """
        write_to_file = False
        self.current_snapshot += 1

        if not self.current_snapshot % 10:
            write_to_file = True
            self.notify_observers(source=__name__,
                                  message="Snapshot #" + str(self.current_snapshot))

        self.simulation.snapshot(write_to_file = write_to_file, 
                                 snapshot_number = self.current_snapshot)
            
    def is_finished(self) -> bool:
        """
        Checks is simulation is finished by checking if maximum number of 
        snashots is reached.
        
        Returns
        -------
            True if simulation is finished; False otherwise.
        """
        if self.current_snapshot < ParametersGeneral.num_snapshots:
            return False
        else:
            return True
            
    def finalize(self):
        """
        Finalizes the simulation and performs all post-simulation tasks
        """
        self.simulation.finalize(snapshot_number=self.current_snapshot)
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
