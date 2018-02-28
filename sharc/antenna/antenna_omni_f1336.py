# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:16:24 2018

@author: Calil
"""

import numpy as np

from sharc.antenna.antenna import Antenna


class AntennaOmniF1336(Antenna):
    """
    Implements an omni antenna as defined in ITU.R F.1336 section 2.2
    (average sidelobe patterns for omnidirectional antennas)
    """
    
    def __init__(self, max_gain: float, down_tilt: float, elevation: float):
        super().__init__()
        # Set attributes
        self.max_gain = max_gain
        self.down_tilt = -1*down_tilt
        self.elevation = elevation
        
        # Fixed parameters
        self.theta_3db = 107.6*(10**(-0.1*max_gain))
        self.k = 0.7
        self.theta_5db = self.theta_3db*np.sqrt(1 - np.log10(self.k + 1)/1.2)
    
    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Theta taken with z axis as reference
        """
        # Get parameters
        theta_vec = np.asarray(kwargs["theta_vec"])
        lo_theta = self.to_local_coord(theta_vec)
        
        theta_e = 90*(lo_theta + self.down_tilt)/(90 + 
                     self.down_tilt*np.sign(lo_theta + self.down_tilt))
        
        theta_abs = np.abs(theta_e)
        gain = np.zeros_like(theta_abs)
        
        # Sections of function
        sec_1 = np.where(np.logical_and(theta_abs >= 0, 
                                        theta_abs < self.theta_3db))
        sec_2 = np.where(np.logical_and(theta_abs >= self.theta_3db, 
                                        theta_abs < self.theta_5db))
        sec_3 = np.where(np.logical_and(theta_abs >= self.theta_5db, 
                                        theta_abs <= 90))
        
        # Calculate for diferent sections
        gain[sec_1] = self.max_gain-12*(theta_e[sec_1]/self.theta_3db)**2
        gain[sec_2] = self.max_gain-15-10*np.log10(self.k+1)
        gain[sec_3] = self.max_gain-15+10*np.log10(
                np.power(theta_abs[sec_3]/self.theta_3db,-1.5)+self.k)
        
        return gain
    
    def to_local_coord(self,theta):
        """
        Receives theta with reference to z axis, converts it to reference in
        x axis and converts to local coordinate system
        """
        lo_theta = np.ravel(np.array([theta + self.elevation]))
        
        lo_theta = 90 - np.ravel(np.mod(np.array([lo_theta]),360))
        
        ofb_theta = np.where(np.logical_or(lo_theta < -90,lo_theta > 90))
        lo_theta[ofb_theta] = np.sign(lo_theta[ofb_theta])*180 - lo_theta[ofb_theta]
        
        return lo_theta
        
    
if __name__ == '__main__':
    # Import
    import matplotlib.pyplot as plt
    
    # Create variables
    max_g = 5
    down_tilt = -10
    elevation = 20
    ant = AntennaOmniF1336(max_g,down_tilt,elevation)
    
    # Gains
    phi_1 = np.linspace(-180, 180, num = 360)
    theta_1 = 90*np.ones_like(phi_1)
    gain_1 = ant.calculate_gain(phi_vec=phi_1,theta_vec=theta_1)
    theta_2 = np.linspace(0,180, num = 360)
    phi_2 = np.zeros_like(theta_2)
    gain_2 = ant.calculate_gain(phi_vec=phi_2,theta_vec=theta_2)
    
    top_y_lim = np.ceil(np.max(gain_1)/10)*10

    fig = plt.figure(figsize=(15,5), facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(121)

    ax1.plot(phi_1,gain_1)
    ax1.set_xlim(-180, 180)
    ax1.grid(True)
    ax1.set_xlabel(r"$\varphi$ [deg]")
    ax1.set_ylabel("Gain [dB]")
    
    ax2 = fig.add_subplot(122, sharey = ax1)
    
    ax2.plot(theta_2,gain_2)
    ax2.set_xlim(0, 180)
    ax2.grid(True)
    ax2.set_xlabel(r"$\theta$ [deg]")
    ax2.set_ylabel("Gain [dB]")
    
    plt.show()
    