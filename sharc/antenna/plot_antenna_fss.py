# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:25:19 2017

@author: edgar
"""

import matplotlib.pyplot as plt
import numpy as np

from sharc.antenna.antenna_fss import AntennaFss
from sharc.parameters.parameters_fss import ParametersFss

if __name__ == '__main__':
    
        # initialize antenna parameters
        param = ParametersFss()
        param.sat_rx_antenna_gain = 51
        param.sat_rx_antenna_pattern = "ITU-R S.672-4"
        param.sat_rx_antenna_l_s = -20    
        param.sat_rx_antenna_3_dB = 0.65
        
        antenna = AntennaFss(param)
        psi = np.linspace(-10,10, num = 360)
        #psi = np.array([0,1])
        antenna.calculate_gain(psi)
    
        #print(antenna.gain)
        
        plt.plot(psi, antenna.gain, "-b")
        plt.ylim((0, 60))
        plt.xlim((-10, 10))
        plt.title("Satellite antenna pattern in the FSS")
        plt.xlabel("angle [degrees]")
        plt.ylabel("Gain [dBi]")
        plt.grid()  
        
        plt.show()