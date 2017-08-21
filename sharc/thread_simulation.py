# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 19:46:20 2017

@author: edgar
"""

from sharc.model import Model

import time
from threading import Thread, Event

class ThreadSimulation(Thread):
    """
    This class extends the Thread class and controls the simulation (start/stop)
    
    Attributes
    ----------
        _stop (Event) : This flag is used to control when simulation is stopped
            by user
        model (Model) : Reference to the Model implementation of MVC
    """
    
    def __init__(self, model: Model):
        Thread.__init__(self)
        self.model = model
        self.stop_flag = Event()
        
    def stop(self):
        """
        This is called by the controller when it receives the stop command by 
        view. This method sets the stop flag that is checked during the 
        simulation. Simulation stops when stop flag is set.
        """
        self.stop_flag.set()
        
    def is_stopped(self) -> bool:
        """
        Checks if stop flag is set.
        
        Returns
        -------
            True if simulation is stopped
        """
        return self.stop_flag.isSet()
        
    def run(self):
        """
        This is overriden from base class and represents the thread's activity.
        """
        start = time.perf_counter()
        
        self.model.initialize()
        while not self.model.is_finished() and not self.is_stopped():
            self.model.snapshot()
        self.model.finalize()
        # calculates simulation time when it finishes and sets the elapsed time
        end = time.perf_counter()
        elapsed_time = time.gmtime(end - start)
        self.model.set_elapsed_time(time.strftime("%H h %M min %S seg", elapsed_time))
        
        
        