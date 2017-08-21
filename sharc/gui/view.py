# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:28:09 2016

@author: edgar
"""

from sharc.controller import Controller
from sharc.support.observer import Observer
from sharc.support.enumerations import Action, State
from sharc.gui.thread_safe_scrolled_text import ThreadSafeScrolledText
from sharc.results import Results

import matplotlib.pyplot as plt

import os
import queue
import logging
import tkinter
import tkinter.filedialog
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
        self.__results_queue = queue.Queue()
        self.__results = None

    def initialize(self):
        """
        Creates all the graphical components
        """
        self.title("SHARC simulator")

        self.__app_icon = tkinter.PhotoImage(file="img/app_icon.gif")
        self.tk.call('wm', 'iconphoto', self._w, self.__app_icon)

        self.__frame = tkinter.Frame(self, bg = '#CCCCCC')
        self.__frame.pack(fill='both', expand='yes')

        self.__scrolledtext = ThreadSafeScrolledText(self.__frame,
            wrap=tkinter.WORD, width=80, height=25, bd=5)
        self.__scrolledtext.grid(column=0, row=1, columnspan=7, sticky='EW')

        self.__start_image = tkinter.PhotoImage(file = "img/start_icon.gif")
        self.__start_button = tkinter.Button(self.__frame, text="START",
            image=self.__start_image, compound=tkinter.LEFT,
            state=tkinter.NORMAL, command=self.__on_start_button_click)
        self.__start_button.bind("<Return>", self.__on_start_button_click)
        self.__start_button.grid(column=1, row=0, sticky='EW')

        self.__stop_image = tkinter.PhotoImage(file="img/stop_icon.gif")
        self.__stop_button = tkinter.Button(self.__frame, text="STOP  ",
            image=self.__stop_image, compound=tkinter.LEFT,
            state=tkinter.DISABLED, command=self.__on_stop_button_click)
        self.__stop_button.bind("<Return>", self.__on_stop_button_click)
        self.__stop_button.grid(column=2, row=0, sticky='EW')

        self.__results_image = tkinter.PhotoImage(file = "img/results_icon.gif")
        self.__results_button = tkinter.Button(self.__frame, text="RESULTS",
            image=self.__results_image, compound=tkinter.LEFT,
            state=tkinter.DISABLED, command=self.__on_results_button_click)
        self.__results_button.bind("<Return>", self.__on_results_button_click)
        self.__results_button.grid(column=3, row=0, sticky='EW')

        self.__clear_image = tkinter.PhotoImage(file="img/clear_icon.gif")
        self.__clear_button = tkinter.Button(self.__frame, text="CLEAR",
            image=self.__clear_image, compound=tkinter.LEFT,
            state=tkinter.NORMAL, command=self.__on_clear_button_click)
        self.__clear_button.grid(column=4, row=0, sticky='EW')

        self.__copy_image = tkinter.PhotoImage(file="img/copy_icon.gif")
        self.__copy_button = tkinter.Button(self.__frame, text="COPY",
            image=self.__copy_image, compound=tkinter.LEFT,
            state=tkinter.NORMAL, command=self.__on_copy_button_click)
        self.__copy_button.grid(column=5, row=0, sticky='EW')

        self.grid_columnconfigure(0, weight=1)
        self.resizable(False, False)
        self.update()
        self.geometry(self.geometry())

        self.__insert_text(__name__, "Ready to run!\n")
        self.__set_state(State.INITIAL)

    def __on_start_button_click(self, *args):
        """
        This method is called when start button is clicked
        """
        default_file = os.path.join(os.getcwd(), "parameters", "parameters.ini")
        default_dir = os.path.join(os.getcwd(), "parameters")
        param_file = tkinter.filedialog.askopenfilename(title = "Select parameters file",
                                                              initialdir = default_dir,
                                                              initialfile = default_file,
                                                              filetypes = (("Simulation parameters", "*.ini"),
                                                                           ("All files", "*.*") ))
        if param_file:
            self.__controller.action(action = Action.START_SIMULATION, 
                                     param_file = param_file)
            

    def __on_stop_button_click(self, *args):
        """
        This method is called when stop button is clicked
        """
        self.__insert_text(__name__, "STOPPED BY USER, FINALIZING SIMULATION")
        self.__controller.action(action = Action.STOP_SIMULATION)
        self.__set_state(State.STOPPING)

    def __on_results_button_click(self, *args):
        """
        This method is called when results button is clicked
        """
        if not self.__results:
            try:
                self.__results = self.__results_queue.get_nowait()
            except queue.Empty:
                pass
        self.__plot_results(self.__results)

    def __on_clear_button_click(self):
        """
        This method is called when clear button is clicked
        """
        self.__scrolledtext.config(state = tkinter.NORMAL)
        self.__scrolledtext.delete(1.0, tkinter.END)
        self.__scrolledtext.config(state = tkinter.DISABLED)

    def __on_copy_button_click(self):
        """
        This method is called when copy button is clicked
        """
        self.clipboard_clear()
        self.__scrolledtext.config(state = tkinter.NORMAL)
        self.clipboard_append(self.__scrolledtext.get(1.0, tkinter.END))
        self.__scrolledtext.config(state = tkinter.DISABLED)
        self.__popup("Copied to clipboard.")

    def __plot_results(self, results: Results):
        file_extension = ".png"
        transparent_figure = False
        
        for plot in results.plot_list:
            plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')
            plt.plot(plot.x, plot.y, color='#990000', linewidth=2)        
            plt.title(plot.title)
            plt.xlabel(plot.x_label)
            plt.ylabel(plot.y_label)
            if not plot.x_lim is None:
                plt.xlim(plot.x_lim)
            if not plot.y_lim is None:
                plt.ylim(plot.y_lim)                
            #plt.grid()
            plt.tight_layout()
            plt.savefig(os.path.join("output", plot.file_name + file_extension), 
                        transparent=transparent_figure)        

        #plt.show()
        self.__popup("Plots successfully created! Check output directory.")

    def __insert_text(self, source: str, text: str):
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
        self.__scrolledtext.write(text + "\n")
        logger = logging.getLogger(source)
        logger.info(text)

    def __set_state(self, state: State):
        """
        Sets the state of the graphical user interface, i.e. enables or
        disables the buttons according to the state of the simulation

        Parameters
        ----------
            state : Enumaration that defines the simulation state
        """
        if state is State.INITIAL:
            self.__start_button.config(state=tkinter.NORMAL)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__results_button.config(state=tkinter.DISABLED)
            self.__start_button.focus_set()
        elif state is State.RUNNING:
            self.__start_button.config(state=tkinter.DISABLED)
            self.__stop_button.config(state=tkinter.NORMAL)
            self.__results_button.config(state=tkinter.NORMAL)
            self.__stop_button.focus_set()
        elif state is State.STOPPING:
            self.__start_button.config(state=tkinter.DISABLED)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__results_button.config(state=tkinter.DISABLED)
        elif state is State.FINISHED or state is State.STOPPED:
            self.__start_button.config(state=tkinter.NORMAL)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__results_button.config(state=tkinter.NORMAL)
        else:
            self.__start_button.config(state=tkinter.NORMAL)
            self.__stop_button.config(state=tkinter.DISABLED)
            self.__results_button.config(state=tkinter.DISABLED)
            self.__start_button.focus_set()

    def set_controller(self, controller: Controller):
        """
        Keeps the reference to the controller

        Parameters
        ----------
            controller : Reference to the controller
        """
        self.__controller = controller

    def notify_observer(self, *args, **kwargs):
        """
        Implements the method from Observer class. See documentation on the
        super class.

        Parameters
        ----------
            state : Enumaration that defines the simulation state
            message : Message that will be displayed on console and on log file
        """
        if "state" in kwargs:
            self.__set_state(kwargs["state"])
            #self.insert_text( __name__, "\n" )
        if "message" in kwargs:
            self.__insert_text(kwargs["source"], kwargs["message"])
        if "results" in kwargs:
            self.__update_results_queue(kwargs["results"])

    def __update_results_queue(self, results: Results):
        self.__results_queue.put(results)

    def __popup(self, message: str):
        """
        Displays a message on a popup window.

        Parameters
        ----------
            message : Message to be displayed on the popup window
        """
        top_level = tkinter.Toplevel(self, bg='#FFFFFF', padx=0, pady=0)
        top_level.title("Message")

        img = tkinter.PhotoImage(file="img/app_icon.gif")
        top_level.tk.call('wm', 'iconphoto', top_level._w, img)

        empty_top_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_top_label.pack()

        label = tkinter.Label(top_level, text=message, height=0, width=50,
                              bg='#FFFFFF')
        label.pack()

        button = tkinter.Button(top_level, text="    OK    ",
                                command=top_level.destroy)
        button.pack()

        empty_botton_label = tkinter.Label(top_level, height=1, bg='#FFFFFF')
        empty_botton_label.pack()

