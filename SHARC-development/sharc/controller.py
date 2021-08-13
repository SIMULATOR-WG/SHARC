# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:32:02 2016

@author: edgar
"""

from sharc.support.enumerations import Action 
from sharc.model import Model

from thread_simulation import ThreadSimulation

class Controller:
    """
    This is the Controller class of the simplified MVC model that is 
    implemented in this application. This class should define application 
    behavior, map user actions to model updates and select view for response.
    """
    
    def __init__(self):
        pass
        
    def set_model(self, model: Model):
        self.model = model

    def get_model(self):
        return self.model
        
    def action(self, *args, **kwargs):
        """
        Receives the user action that is captured by view and maps it to the 
        appropriate action. Currently, the only supported actions are the ones
        that start and stop simulation. 
        
        Parameters
        ----------
            Action: this non-keyworded argument indicates the action to be taken
        """
        action = kwargs["action"]
        
        if action is Action.START_SIMULATION:
            self.model.set_param_file(kwargs["param_file"])
            self.simulation_thread = ThreadSimulation(self.model)
            self.simulation_thread.start()
        if action is Action.START_SIMULATION_SINGLE_THREAD:
            self.model.set_param_file(kwargs["param_file"])
            self.simulation_thread = ThreadSimulation(self.model)
            # call run method directly, without starting a new thread
            self.simulation_thread.run()        
        if action is Action.STOP_SIMULATION:
            if self.simulation_thread.is_alive():
                self.simulation_thread.stop()
                