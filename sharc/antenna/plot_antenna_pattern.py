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
    
    def plot_element_pattern(self,antenna: AntennaBeamformingImt, sta_type: str, antenna_type: str, plot_type: str):
        
        phi_escan = 0
        theta_tilt = 90
        
        # Plot horizontal pattern
        phi = np.linspace(-180,180, num = 360)
        theta = theta_tilt*np.ones(np.size(phi))

        if plot_type == "ELEMENT":
            gain = antenna.element.element_pattern(phi, theta)
        elif plot_type == "ARRAY":
            antenna.add_beam(phi_escan,theta_tilt)
            gain = antenna.calculate_gain(phi,theta,np.zeros_like(phi, dtype=int))

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
            
        ax1.set_xlim(-180, 180)

        # Plot vertical pattern
        theta = np.linspace(0,180, num = 360)
        phi = (0 + phi_escan)*np.ones(np.size(theta))

        if plot_type == "ELEMENT":
            gain = antenna.element.element_pattern(phi, theta)
        elif plot_type == "ARRAY":
            gain = antenna.calculate_gain(phi,theta,np.zeros_like(phi, dtype=int))

        ax2 = fig.add_subplot(122, sharey = ax1)

        ax2.plot(theta,gain)
        ax2.grid(True)
        ax2.set_xlabel(r"$\theta$ [deg]")
        ax2.set_ylabel("Gain [dB]")
        
        if plot_type == "ELEMENT":
            ax2.set_title("Element Vertical Radiation Pattern")
        elif plot_type == "ARRAY":
            ax2.set_title("Array Vertical Radiation Pattern")
        
        ax2.set_xlim(0, 180)
        top_y_lim = np.ceil(np.max(gain)/10)*10
        ax2.set_ylim(top_y_lim - 50,top_y_lim)
        
        if sta_type == "BS":
            file_name = self.figs_dir + "bs_"
        elif sta_type == "UE":
            file_name = self.figs_dir + "ue_"
            
        if antenna_type == "TX":
            file_name = file_name + "tx_"
        elif antenna_type == "RX":
            file_name = file_name + "rx_"
            
        if plot_type == "ELEMENT":
            file_name = file_name + "element_pattern.png"
        elif plot_type == "ARRAY":
            file_name = file_name + "array_pattern.png"
        
        plt.savefig(file_name)
        
        
if __name__ == '__main__':
    
    figs_dir = "figs/"

    param = ParametersAntennaImt()
    plot = PlotAntennaPattern(figs_dir)

    # Plot BS TX radiation patterns
    par = param.get_antenna_parameters("BS","TX")
    bs_array = AntennaBeamformingImt(par,0,0)
    plot.plot_element_pattern(bs_array,"BS","TX","ELEMENT")
    plot.plot_element_pattern(bs_array,"BS","TX","ARRAY")
    
    # Plot UE TX radiation patterns
    par = param.get_antenna_parameters("UE","TX")
    ue_array = AntennaBeamformingImt(par,0,0)
    plot.plot_element_pattern(ue_array,"UE","TX","ELEMENT")
    plot.plot_element_pattern(ue_array,"UE","TX","ARRAY")
    
    # Plot BS RX radiation patterns
    par = param.get_antenna_parameters("BS","RX")
    bs_array = AntennaBeamformingImt(par,0,0)
    plot.plot_element_pattern(bs_array,"BS","RX","ELEMENT")
    plot.plot_element_pattern(bs_array,"BS","RX","ARRAY")
    
    # Plot UE RX radiation patterns
    par = param.get_antenna_parameters("UE","RX")
    ue_array = AntennaBeamformingImt(par,0,0)
    plot.plot_element_pattern(ue_array,"UE","RX","ELEMENT")
    plot.plot_element_pattern(ue_array,"UE","RX","ARRAY")
    
    print('END')
    