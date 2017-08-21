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

from sharc.antenna.antenna_fss_9_1m import AntennaFss_9_1

class PlotAntennaPatternFSS_9_1(object):
    """"
    
    
    """
    
    def __init__(self,gain_dir):
        self.gain_dir = gain_dir
        
    def plot_pattern(self,antenna: AntennaFss_9_1):
                
        # Plot radiation pattern for theta = 90 degrees
        phi = np.arange (0,10.1,0.1)
        gain = np.array(antenna.get_gain())

        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(121)
        ax1.plot(phi,gain)
        ax1.grid(True)
        ax1.set_xlabel(r"$\varphi$ [deg]")
        ax1.set_ylabel("Gain [dBi]")
        ax1.set_title("Radiation Pattern")     
        ax1.set_xlim(0, 10)

        # Plot radiation pattern for theta = 45 degrees
        
        file_name = self.gain_dir + "FSS_9_1.png"

        self.save_file(file_name)
        
    def save_file(self,file_name: str):
        plt.savefig(file_name)
        
if __name__ == '__main__':
    
    gain_dir = "/Users/carlosrodriguez/Desktop/"

    plot = PlotAntennaPatternFSS_9_1(gain_dir)
    antenna_array = AntennaFss_9_1(62)
    plot.plot_pattern(antenna_array)

