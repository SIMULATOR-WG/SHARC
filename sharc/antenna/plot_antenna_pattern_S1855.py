#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class and script for plotting antenna pattern for Earth Station antenna for FSS
service. 

Created on Mon Jun  5 10:31:24 2017

@author: carlosrodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

from sharc.antenna.antenna_s1855 import AntennaS1855

class PlotAntennaPatternS1855(object):
    """"
    Implements plot of the antenna radiation pattern for Earth Station according
    Reccomendation  ITU-R S.1855 (01/2010)
    "Alternative reference radiation pattern for earth station antennas used with
    satellites in the geostationary-satellite orbit for use in coordination and/or
    interference assesment in the frequency range from 2 to 31 GHz"
    
    """
    
    def __init__(self,gain_dir):
        self.gain_dir = gain_dir
        
    def plot_pattern(self,antenna: AntennaS1855):
                
        # Plot radiation pattern for theta = 90 degrees
        phi = np.linspace(0,180, num = 180)
        theta = np.array([90.0])

        gain = np.array(antenna.get_gain(phi,theta))

        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(121)
        ax1.plot(phi,gain)
        ax1.grid(True)
        ax1.set_xlabel(r"$\varphi$ [deg]")
        ax1.set_ylabel("Gain [dBi]")
        ax1.set_title("Radiation Pattern $\theta$=90 degrees")     
        ax1.set_xlim(0, 180)

        # Plot radiation pattern for theta = 45 degrees
        phi = np.linspace(0,180, num = 180)
        theta = np.array([45.0])

        gain = np.array(antenna.get_gain(phi,theta))

        ax2 = fig.add_subplot(122, sharey = ax1)
        ax2.plot(phi,gain)
        ax2.grid(True)
        ax2.set_xlabel(r"$\varphi$ [deg]")
        ax2.set_ylabel("Gain [dBi]")
        
        ax2.set_title("Radiation Pattern $\theta$=45 degrees")

        
        ax2.set_xlim(0, 180)
        top_y_lim = np.ceil(np.max(gain)/10)*10
        ax2.set_ylim(top_y_lim - 100,top_y_lim)
        
        file_name = self.gain_dir + "S1855.png"

        self.save_file(file_name)
        
    def save_file(self,file_name: str):
        plt.savefig(file_name)
        
if __name__ == '__main__':
    
    gain_dir = "/Users/carlosrodriguez/Desktop/"

    plot = PlotAntennaPatternS1855(gain_dir)
    antenna_array = AntennaS1855(9.1, 24250, 62)
    plot.plot_pattern(antenna_array)

