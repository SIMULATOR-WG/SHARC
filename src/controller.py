# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:32:02 2016

@author: edgar
"""

from support.enumerations import Action 
from simulation_thread import SimulationThread

class Controller:
    """
    This is the Controller class of the simplified MVC model that is 
    implemented in this application. This class should define application 
    behavior, map user actions to model updates and select view for response.
    """
    
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.model.add_observer(view)
        self.view.add_controller(self)
        
    def action(self, *args):
        """
        Receives the user action that is captured by view and maps it to the 
        appropriate action. Currently, the only supported actions are the ones
        that start and stop simulation. 
        
        Parameters
        ----------
        Action
            This non-keyworded argument indicates the action to be taken
            
        """
        if Action.START_SIMULATION in args:
            self.simulation_thread = SimulationThread(self.model)
            self.simulation_thread.start()
            #self.simulation_thread.run()
        if Action.STOP_SIMULATION in args:
            if self.simulation_thread.is_alive():
                self.simulation_thread.stop()
                