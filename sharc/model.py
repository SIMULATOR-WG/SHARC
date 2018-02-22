# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 17:03:51 2016

@author: edgar
"""

from sharc.support.observable import Observable
from sharc.support.observer import Observer
from sharc.support.enumerations import State
from sharc.simulation_downlink import SimulationDownlink
from sharc.simulation_uplink import SimulationUplink
from sharc.parameters.parameters import Parameters

import random
import sys

class Model(Observable):
    """
    Implements the Observable interface. It has a reference to the simulation
    object and controls the simulation flow (init/step/finilize).
    """

    def __init__(self):
        super(Model, self).__init__()
        self.simulation = None
        self.parameters = None
        self.param_file = None

    def add_observer(self, observer: Observer):
        Observable.add_observer(self, observer)

    def set_param_file(self, param_file):
        self.param_file = param_file
        self.notify_observers(source = __name__,
                              message = "Loading file:\n" + self.param_file)

    def initialize(self):
        """
        Initializes the simulation and performs all pre-simulation tasks, such
        as loading parameters.
        """
        self.parameters = Parameters()
        self.parameters.set_file_name(self.param_file)
        self.parameters.read_params()

        if self.parameters.general.imt_link == "DOWNLINK":
            self.simulation = SimulationDownlink(self.parameters, self.param_file)
        else:
            self.simulation = SimulationUplink(self.parameters, self.param_file)
        self.simulation.add_observer_list(self.observers)

        description = self.get_description()

        self.notify_observers(source=__name__,
                              message=description + "\nSimulation is running...",
                              state=State.RUNNING )
        self.current_snapshot = 0

        self.simulation.initialize()

        random.seed( self.parameters.general.seed )

        self.secondary_seeds = [None] * self.parameters.general.num_snapshots

        max_seed = 2**32 - 1

        for index in range(self.parameters.general.num_snapshots):
            self.secondary_seeds[index] = random.randint(1, max_seed)

    def get_description(self) -> str:
        param_system = self.simulation.param_system

        description = "\nIMT:\n" \
                            + "\tinterfered with: {:s}\n".format(str(self.parameters.imt.interfered_with)) \
                            + "\tdirection: {:s}\n".format(self.parameters.general.imt_link) \
                            + "\tfrequency: {:.3f} GHz\n".format(self.parameters.imt.frequency*1e-3) \
                            + "\tbandwidth: {:.0f} MHz\n".format(self.parameters.imt.bandwidth) \
                            + "\ttopology: {:s}\n".format(self.parameters.imt.topology) \
                            + "\tpath loss model: {:s}\n".format(self.parameters.imt.channel_model)  \
                    + "{:s}:\n".format(self.parameters.general.system) \
                            + "\tfrequency: {:.3f} GHz\n".format(param_system.frequency*1e-3) \
                            + "\tbandwidth: {:.0f} MHz\n".format(param_system.bandwidth) \
                            + "\tpath loss model: {:s}\n".format(param_system.channel_model) \
                            + "\tantenna pattern: {:s}\n".format(param_system.antenna_pattern)

        return description

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

        self.simulation.snapshot(write_to_file=write_to_file,
                                 snapshot_number=self.current_snapshot,
                                 seed = self.secondary_seeds[self.current_snapshot - 1])

    def is_finished(self) -> bool:
        """
        Checks is simulation is finished by checking if maximum number of
        snashots is reached.

        Returns
        -------
            True if simulation is finished; False otherwise.
        """
        if self.current_snapshot < self.parameters.general.num_snapshots:
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
