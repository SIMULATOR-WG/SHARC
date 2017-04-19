# -*- coding: utf-8 -*-
"""
Class and script for plotting antenna patterns. Script in the end of file.
Simply run this file to plot patterns.

Created on Sun Apr 16 17:49:16 2017

@author: Calil
"""

import numpy as np
import matplotlib.pyplot as plt

from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt

class PlotAntennaPattern(object):
    
    def __init__(self,figs_dir):
        self.figs_dir = figs_dir
    
    def plot_element_pattern(self,antenna: AntennaBeamformingImt, sta_type: str, plot_type: str):
        
        phi_escan = 180
        theta_tilt = 0
        
        # Plot horizontal pattern
        phi = np.linspace(-180,180, num = 360)
        theta = 90 + theta_tilt

        if plot_type == "ELEMENT":
            gain = np.array([antenna.element_pattern(p,theta) for p in phi])
        elif plot_type == "ARRAY":
            gain = np.array([antenna.beam_gain(p,theta,phi_escan,theta_tilt) for p in phi])

        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(121)

        ax1.plot(phi,gain)
        ax1.grid(True)
        ax1.set_xlabel(r"$\varphi$ [deg]")
        ax1.set_ylabel("Gain [dB]")
        
        if plot_type == "ELEMENT":
            ax1.set_title("Element Horizontal Radiation Pattern")
        elif plot_type == "ARRAY":
            ax1.set_title("Array Horizontal Radiation Pattern")
            
#        ax1.set_xlim(-180, 180)

        # Plot vertical pattern
        phi = 0 + phi_escan
        theta = np.linspace(0,180, num = 360)

        if plot_type == "ELEMENT":
            gain = np.array([antenna.element_pattern(phi,t) for t in theta])
        elif plot_type == "ARRAY":
            gain = np.array([antenna.beam_gain(phi,t,phi_escan,theta_tilt) for t in theta])

        ax2 = fig.add_subplot(122, sharey = ax1)

        ax2.plot(theta,gain)
        ax2.grid(True)
        ax2.set_xlabel(r"$\theta$ [deg]")
        ax2.set_ylabel("Gain [dB]")
        
        if plot_type == "ELEMENT":
            ax2.set_title("Element Vertical Radiation Pattern")
        elif plot_type == "ARRAY":
            ax2.set_title("Array Vertical Radiation Pattern")
        
#        ax2.set_xlim(0, 180)
#        top_y_lim = np.ceil(np.max(gain)/10)*10
#        ax2.set_ylim(top_y_lim - 50,top_y_lim)
        
        if sta_type == "BS":
            file_name = self.figs_dir + "bs_"
        elif sta_type == "UE":
            file_name = self.figs_dir + "ue_"
            
        if plot_type == "ELEMENT":
            file_name = file_name + "element_pattern.png"
        elif plot_type == "ARRAY":
            file_name = file_name + "array_pattern.png"
        
        self.save_file(file_name)
        
    def save_file(self,file_name: str):
        plt.savefig(file_name)
        
if __name__ == '__main__':
    
    figs_dir = ""

    param = ParametersAntennaImt()
    plot = PlotAntennaPattern(figs_dir)

    # Plot BS radiation patterns
    bs_array = AntennaBeamformingImt(param,"BS")
    plot.plot_element_pattern(bs_array,"BS","ELEMENT")
    plot.plot_element_pattern(bs_array,"BS","ARRAY")
    
    # Plot UE radiation patterns
    ue_array = AntennaBeamformingImt(param,"UE")
    plot.plot_element_pattern(ue_array,"UE","ELEMENT")
    plot.plot_element_pattern(ue_array,"UE","ARRAY")
    