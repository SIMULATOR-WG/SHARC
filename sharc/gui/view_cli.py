# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:18:43 2017

@author: edgar
"""
from sharc.controller import Controller
from sharc.support.observer import Observer
from sharc.support.enumerations import Action

import logging


class ViewCli(Observer):
    """
    Implements the command line interface. This is a subclass of Observer and
    has to implement the notify_observer method.
    """

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent

        
    def initialize(self, param_file):
        self.controller.action(action = Action.START_SIMULATION_SINGLE_THREAD, 
                                     param_file = param_file)
        
        
    def set_controller(self, controller: Controller):
        """
        Keeps the reference to the controller

        Parameters
        ----------
            controller : Reference to the controller
        """
        self.controller = controller

        
    def notify_observer(self, *args, **kwargs):
        """
        Implements the method from Observer class. See documentation on the
        super class.

        Parameters
        ----------
            state : Enumaration that defines the simulation state
            message : Message that will be displayed on console and on log file
        """
        if "message" in kwargs:
            self.insert_text(kwargs["source"], kwargs["message"])

            
    def insert_text(self, source: str, text: str):
        """
        This method can be called to display a message on the console. The same
        message will be written to the log file and the file will also include
        the name of the class that called the method. By default, all logging
        messages will be written in INFO level.

        Parameters
        ----------
            source : Class name that called the method
            text : Message that will be displayed on console and on log file
        """
        logger = logging.getLogger(source)
        logger.info(text)
        