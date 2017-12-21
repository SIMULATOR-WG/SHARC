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
    
    def __init__(self, gain, azimuth: float, elevation: float):
        """
        Constructs an AntennaM1466 object.
        Does not receive angles in local coordinate system.
        Elevation taken wrt x-y plane.
        
        Parameters
        ---------
            azimuth (float): antenna's physical azimuth inclination
            elevation (float): antenna's physical elevation inclination
                wrt x-y plane
        """
        super().__init__()
        self.peak_gain = gain
        self.azimuth = azimuth
        self.elevation = elevation

        
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
        
        phi_l, theta_l = self.to_local_coord(phi, theta)
        
        gain = self.peak_gain + self.get_gain_az(phi_l) + self.get_gain_elev(theta_l)
        
        return gain
        
        
    def get_gain_az(self, phi: np.array) -> np.array:
        """
        Returns the antenna gain in the azimuth plane for the given direction.
        """
        pa = np.abs(phi)
        gain = np.empty(pa.shape)
        gain[pa < 0.5] = 0
        gain[(pa >= 0.5) & (pa < 2)] = -10
        gain[(pa >= 2) & (pa < 5)] = -20
        gain[(pa >= 5) & (pa < 40)] = -27.5
        gain[(pa >= 40)] = -35
        return gain
    
    
    def get_gain_elev(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the elevation plane for the given direction.
        """
        gain = np.empty(theta.shape)
        gain[theta < -60] = -35
        gain[(theta >= -60) & (theta < -30)] = -27.5
        gain[(theta >= -30) & (theta < 5)] = 0
        gain[(theta >= 5) & (theta < 10)] = -20
        gain[theta >= 10] = -35
        return gain
        
        
    def to_local_coord(self, phi: np.array, theta: np.array) -> tuple:
        """
        Converts the azimuth and elevation angles to the local coordinate system,
        taking into account the physical orientation of the antenna.
        """
        
        phi_l = phi - self.azimuth
        theta_l = 90 - theta - self.elevation
        
        return phi_l, theta_l

        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    azimuth = np.linspace(-90, 90, num = 100000)
    elevation = np.linspace(-80, 20, num = 100000)
    
    gain = 30
    azimuth_l = 0
    elevation_l = 0
    
    antenna = AntennaM1466(gain, azimuth_l, elevation_l)
    
    gain_az = gain + antenna.get_gain_az(azimuth)
    gain_elev = gain + antenna.get_gain_elev(elevation)

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