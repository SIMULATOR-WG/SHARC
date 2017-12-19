# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 16:15:18 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna

import numpy as np

class AntennaM1466(Antenna):
    """
    Implements the antenna pattern for radionavigation radars according to 
    Rec. ITU-R M.1466-1
    """
    
    def __init__(self):
        super().__init__()
        self.peak_gain = 30

        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the gain in the given direction.
        
        Parameters
        ----------
            phi_vec (np.array): azimuth angles [degrees]
            theta_vec (np.array): elevation angles [degrees]
            
        Returns
        -------
            gain (np.array): gain corresponding to each of the given directions
        """
        phi = np.asarray(kwargs["phi_vec"])
        theta = np.asarray(kwargs["theta_vec"])        
        
        gain = self.peak_gain + self.get_gain_az(phi) + self.get_gain_elev(theta)
        
        return gain
        
        
    def get_gain_az(self, phi: np.array) -> np.array:
        pa = np.abs(phi)
        gain = np.empty(pa.shape)
        gain[pa < 0.5] = 0
        gain[(pa >= 0.5) & (pa < 2)] = -10
        gain[(pa >= 2) & (pa < 5)] = -20
        gain[(pa >= 5) & (pa < 40)] = -27.5
        gain[(pa >= 40)] = -35
        return gain
    
    
    def get_gain_elev(self, theta: np.array) -> np.array:
        gain = np.empty(theta.shape)
        gain[theta < -60] = -35
        gain[(theta >= -60) & (theta < -30)] = -27.5
        gain[(theta >= -30) & (theta < 5)] = 0
        gain[(theta >= 5) & (theta < 10)] = -20
        gain[theta >= 10] = -35
        return gain
        
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    azimuth = np.linspace(-90, 90, num = 100000)
    elevation = np.linspace(-80, 20, num = 100000)
    
    antenna = AntennaM1466()
    
    gain_az = antenna.calculate_gain(phi_vec = azimuth, theta_vec = 0)
    gain_elev = antenna.calculate_gain(phi_vec = 0, theta_vec = elevation)

    fig = plt.figure(figsize=(15,5), facecolor='w', edgecolor='k')
    #plt.title("Rec. ITU-R M.1466 antenna pattern")
    ax1 = fig.add_subplot(121)
    ax1.plot(azimuth, gain_az)
    ax1.grid(True)
    ax1.set_xlabel("Azimuth angle [deg]")
    ax1.set_ylabel("Antenna gain [dBi]")
    ax1.set_xlim([-90, 90])
    ax1.set_ylim([-10, 35])
    
    ax2 = fig.add_subplot(122)
    ax2.plot(elevation, gain_elev)
    ax2.grid(True)
    ax2.set_xlabel("Elevation angle [deg]")
    ax2.set_ylabel("Antenna gain [dBi]")
    ax2.set_xlim([-80, 20])
    ax2.set_ylim([-10, 35])

    plt.show()        