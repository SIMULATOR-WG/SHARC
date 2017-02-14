# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:28:09 2016

@author: edgar
"""

from controller import Controller
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

    def __init__(self, parent=None):
        super(View, self).__init__()
        self.parent = parent
        self.initialize()
        
    def initialize(self):
        """
        Creates all the graphical components
        """
        self.title("SHARC simulator")
        
        img = tkinter.PhotoImage(file="img/app_icon.png")
        self.tk.call('wm', 'iconphoto', self._w, img)        
        
        self.__frame = tkinter.Frame(self, bg = '#CCCCCC')
        self.__frame.pack(fill='both', expand='yes')

        self.__scrolledtext = tkinter.scrolledtext.ScrolledText(self.__frame, 
            wrap=tkinter.WORD, width=80, height=25, bd=5)
        self.__scrolledtext.grid(column=0, row=1, columnspan=2, sticky='EW')

        self.__start_image = tkinter.PhotoImage(file = "img/start_icon.png")
        self.__start_button = tkinter.Button(self.__frame, text=" START",  
            image=self.__start_image, compound=tkinter.LEFT, 
            state=tkinter.NORMAL, command=self.__on_start_button_click)
        self.__start_button.bind("<Return>", self.__on_start_button_click)
        self.__start_button.grid(column=0, row=0, sticky='E')
        
        self.__stop_image = tkinter.PhotoImage(file="img/stop_icon.png")
        self.__stop_button = tkinter.Button(self.__frame, text=" STOP  ",  
            image=self.__stop_image, compound=tkinter.LEFT, 
            state=tkinter.DISABLED, command=self.__on_stop_button_click)
        self.__stop_button.bind("<Return>", self.__on_stop_button_click)
        self.__stop_button.grid(column=1, row=0, sticky='W')
        
        self.__copy_image = tkinter.PhotoImage(file="img/copy_icon.png")
        self.__copy_button = tkinter.Button(self.__frame, text=" COPY  ",  
            image=self.__copy_image, compound=tkinter.LEFT, 
            state=tkinter.NORMAL, command=self.__on_copy_button_click)
        self.__copy_button.grid(column=1, row=0, sticky='E')        
        
        self.grid_columnconfigure(0, weight=1)
        self.resizable(False, False)
        self.update()
        self.geometry(self.geometry())    
        
        self.__insert_text(__name__, "Ready to run!")
        self.__set_state(State.INITIAL)
        
    def __on_start_button_click(self, *args):
        """
        This method is called when start button is clicked
        """
        self.__controller.action(Action.START_SIMULATION)
        
    def __on_stop_button_click(self, *args):
        """
        This method is called when stop button is clicked
        """
        self.__insert_text(__name__, "STOPPED BY USER, FINALIZING SIMULATION")
        self.__controller.action(Action.STOP_SIMULATION)
        self.__set_state(State.STOPPING)
        
    def __on_copy_button_click(self):
        """
        This method is called when copy button is clicked
        """        
        self.clipboard_clear()
        self.clipboard_append(self.__scrolledtext.get(1.0, tkinter.END))  
        self.__popup("Copied to clipboard.")
        
    def __insert_text(self, source: str, text: str):
        """
        This method can be called to display a message on the console. The same
        message will be written to the log file and the file will also include 
        the name of the class that called the method. By default, all logging
        messages will be written in INFO level.
        
        Args:
            source : Class name that called the method
            text : Message that will be displayed on console and on log file
        """
        self.__scrolledtext.insert(tkinter.INSERT, text + "\n")
        self.__scrolledtext.see(tkinter.END)
        logger = logging.getLogger(source)
        logger.info(text)
        
    def __set_state(self, state: State):
        """
        Sets the state of the graphical user interface, i.e. enables or
        disables the buttons according to the state of the simulation
        
        Args:
            state : Enumaration that defines the simulation state
        """
        if state == State.INITIAL or state == State.FINISHED or state == State.STOPPED:
            self.__start_button.config(state=tkinter.NORMAL)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__start_button.focus_set()
        elif state == State.RUNNING:
            self.__start_button.config(state=tkinter.DISABLED)
            self.__stop_button.config(state=tkinter.NORMAL)
            self.__stop_button.focus_set()
        elif state == State.STOPPING:
            self.__start_button.config(state=tkinter.DISABLED)
            self.__stop_button.config(state=tkinter.DISABLED)
        else:
            self.__start_button.config(state=tkinter.NORMAL)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__start_button.focus_set()
            
    def set_controller(self, controller: Controller):
        """
        Keeps the reference to the controller
        
        Args:
            controller : Reference to the controller
        """
        self.__controller = controller
        
    def notify_observer(self, *args, **kwargs):
        """
        Implements the method from Observer class. See documentation on the 
        super class.
        
        Args:
            state : Enumaration that defines the simulation state
            message : Message that will be displayed on console and on log file
        """
        if "state" in kwargs:
            self.__set_state(kwargs["state"])
            #self.insert_text( __name__, "\n" )
        if "message" in kwargs:
            self.__insert_text(kwargs["source"], kwargs["message"])
        
    def __popup(self, message: str):
        """
        Displays a message on a popup window.
        
        Args:
            message : Message to be displayed on the popup window
        """
        top_level = tkinter.Toplevel(self, bg='#FFFFFF', padx=0, pady=0)
        top_level.title("Message")
        
        img = tkinter.PhotoImage(file="img/app_icon.png")
        top_level.tk.call('wm', 'iconphoto', top_level._w, img)    
        
        empty_top_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_top_label.pack()
        
        label = tkinter.Label(top_level, text=message, height=0, width=35, 
                              bg='#FFFFFF')
        label.pack()
        
        button = tkinter.Button(top_level, text="    OK    ", 
                                command=top_level.destroy)
        button.pack()
        
        empty_botton_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_botton_label.pack()
       
