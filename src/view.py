# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:28:09 2016

@author: edgar
"""

from support.observer import Observer
from support.enumerations import Action, State

import logging
import tkinter
import tkinter.scrolledtext

class View(tkinter.Tk, Observer):
    """
    Implements the graphical user interface. This is a subclass of Observer and
    has to implement the notify_observer method.

    """    

    def __init__(self, parent):
        tkinter.Tk.__init__(self, parent)
        self.parent = parent
        self.initialize()
        
    def initialize(self):
        """
        Creates all the graphical components
        """
        self.title("Simulator")
        
        self.__frame = tkinter.Frame(self, bg = '#CCCCCC')
        self.__frame.pack(fill='both', expand='yes')

        self.__scrolledtext = tkinter.scrolledtext.ScrolledText(self.frame, 
            wrap=tkinter.WORD, width=80, height=25, bd=5)
        self.scrolledtext.grid(column=0, row=1, columnspan=2, sticky='EW')

        self.start_image = tkinter.PhotoImage(file = "img/start_icon.png")
        self.start_button = tkinter.Button(self.frame, text=" START",  
            image=self.start_image, compound=tkinter.LEFT, 
            state=tkinter.NORMAL, command=self.on_start_button_click)
        self.start_button.bind("<Return>", self.on_start_button_click)
        self.start_button.grid(column=0, row=0, sticky='E')
        
        self.stop_image = tkinter.PhotoImage(file="img/stop_icon.png")
        self.stop_button = tkinter.Button(self.frame, text=" STOP  ",  
            image=self.stop_image, compound=tkinter.LEFT, 
            state=tkinter.DISABLED, command=self.on_stop_button_click)
        self.stop_button.bind("<Return>", self.on_stop_button_click)
        self.stop_button.grid(column=1, row=0, sticky='W')
        
        self.copy_image = tkinter.PhotoImage(file="img/copy_icon.png")
        self.copy_button = tkinter.Button(self.frame, text=" COPY  ",  
            image=self.copy_image, compound=tkinter.LEFT, 
            state=tkinter.NORMAL, command=self.on_copy_button_click)
        self.copy_button.grid(column=1, row=0, sticky='E')        
        
        self.grid_columnconfigure(0, weight=1)
        self.resizable(False, False)
        self.update()
        self.geometry(self.geometry())    
        
        self.insert_text(__name__, "Ready to run!")
        self.set_state(State.INITIAL)
        
    def on_start_button_click(self, *args):
        """
        This method is called when start button is clicked
        """
        self.controller.action(Action.START_SIMULATION)
        
    def on_stop_button_click(self, *args):
        """
        This method is called when stop button is clicked
        """
        self.insert_text(__name__, "STOPPED BY USER, FINALIZING SIMULATION")
        self.controller.action(Action.STOP_SIMULATION)
        self.set_state(State.STOPPING)
        
    def on_copy_button_click(self):
        """
        This method is called when copy button is clicked
        """        
        self.clipboard_clear()
        self.clipboard_append(self.scrolledtext.get(1.0, tkinter.END))  
        self.popup("Copied to clipboard.")
        
    def insert_text(self, source, text):
        """
        This method can be called to display a message on the console. The same
        message will be written to the log file and the file will also include 
        the name of the class that called the method. By default, all logging
        messages will be written in INFO level.
        
        Parameters
        ----------
            source : string
                Class name that called the method
            text : string
                Message that will be displayed on console and on log file
        """
        self.scrolledtext.insert(tkinter.INSERT, text + "\n")
        self.scrolledtext.see(tkinter.END)
        logger = logging.getLogger(source)
        logger.info(text)
        
    def set_state(self, state):
        """
        Sets the state of the graphical user interface, i.e. enables or
        disables the buttons according to the state of the simulation
        
        Parameter
        ---------
            state : State
                Enumaration that defines the simulation state
        """
        if state == State.INITIAL or state == State.FINISHED or state == State.STOPPED:
            self.start_button.config(state=tkinter.NORMAL)
            self.stop_button.config(state=tkinter.DISABLED)
            self.start_button.focus_set()
        elif state == State.RUNNING:
            self.start_button.config(state=tkinter.DISABLED)
            self.stop_button.config(state=tkinter.NORMAL)
            self.stop_button.focus_set()
        elif state == State.STOPPING:
            self.start_button.config(state=tkinter.DISABLED)
            self.stop_button.config(state=tkinter.DISABLED)
        else:
            self.start_button.config(state=tkinter.NORMAL)
            self.stop_button.config(state=tkinter.DISABLED)
            self.start_button.focus_set()
            
    def add_controller(self, controller):
        """
        Keeps the reference to the controller
        
        Parameter
        ---------
            controller : Controller
                Reference to the controller
        """
        self.controller = controller
        
    def notify_observer(self, *args, **kwargs):
        """
        Implements the method from Observer class. See documentation on the 
        super class.
        
        """
        if "state" in kwargs:
            self.set_state(kwargs["state"])
            #self.insert_text( __name__, "\n" )
        if "message" in kwargs:
            self.insert_text(kwargs["source"], kwargs["message"])
        
    def popup(self, message):
        """
        Displays a message on a popup window.
        
        Parameter
        ---------
            message : string
                Message to be displayed on the popup window
        
        """
        top_level = tkinter.Toplevel(self, bg='#FFFFFF', padx=0, pady=0)
        top_level.title("Message")
        
        empty_top_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_top_label.pack()
        
        label = tkinter.Label(top_level, text=message, height=0, width=35, 
                              bg='#FFFFFF')
        label.pack()
        
        button = tkinter.Button(top_level, text="    OK    ", command=top_level.destroy)
        button.pack()
        
        empty_botton_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_botton_label.pack()
       
