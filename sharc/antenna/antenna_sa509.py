<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 15:09:46 2017

@author: carlosrodriguez
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs
=======
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:34:11 2017

@author: Calil
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_ras import ParametersRas
>>>>>>> dev-ras

import numpy as np

class AntennaSA509(Antenna):
    """
<<<<<<< HEAD
    Implements the Earth station antenna pattern in the SRS service
    according to Recommendation ITU-R SA.509-3
    
    Implementation for multiple entry interference item 1.2 of Rec.
    Recommendation for D/lamnda >=100 always
    """
    
    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 3e8 / ( param.frequency * 1e6 )

        
        self.phi_0 = 20 * 1.732 * lmbda/param.diameter 
        self.phi_1 = self.phi_0 * 2.582 
        self.phi_2 = np.power(10,(49-self.peak_gain)/25)
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros(phi.shape)
                
        idx_0 = np.where(phi < self.phi_1)[0]
        gain[idx_0] = self.peak_gain - 3 * (phi[idx_0]/self.phi_0) * (phi[idx_0]/self.phi_0)
        
        idx_1 = np.where((self.phi_1 <= phi) & (phi < self.phi_2))[0]
        gain[idx_1] = self.peak_gain - 20

        idx_2 = np.where((self.phi_2 <= phi) & (phi < 48))[0]
        gain[idx_2] = 29 - 25 * np.log10(phi[idx_2])

        idx_3 = np.where((48 <= phi) & (phi <= 80))[0]
        gain[idx_3] = -13
        
        idx_4 = np.where((80 <= phi) & (phi <= 120))[0]
        gain[idx_4] = -8
        
        idx_5 = np.where((120 <= phi) & (phi <= 180))[0]
        gain[idx_5] = -13
        
        return gain
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 100, num = 100000)
    
    # initialize antenna parameters
    param27 = ParametersFssEs()
    
 
    param27.antenna_pattern = "ITU-R S.509-3"
    param27.frequency = 23800
    param27.antenna_gain = 50
    param27.diameter = 9.6
    antenna27 = AntennaSA509(param27)
    

    gain27 = antenna27.calculate_gain(phi_vec=phi)
    

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.semilogx(phi, gain27 - param27.antenna_gain, "-b", label = "$f = 23.8$ $GHz,$ $D = 9.6$ $m$")

    plt.title("ITU-R SA.509 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-80, 10))
    
    #ax = plt.gca()
    #ax.set_yticks([-30, -20, -10, 0])
    #ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()        
=======
    Implements the antenna pattern for the Radio Astronomy Service 
    according to recommendation ITU-R SA.509-3.
    """
    
    def __init__(self, param: ParametersRas):
        super().__init__()
        # Set basic attributes
        self.diameter = param.diameter
        self.efficiency = param.antenna_efficiency
        self.wavelength = param.SPEED_OF_LIGHT/(param.frequency*1e6)
        
        # Effective area
        self.effective_area = self.efficiency*(np.pi*self.diameter**2)/4
        
        # Diagram parameters
        self.g_0 = 10*np.log10(self.efficiency*\
                               (np.pi*self.diameter/self.wavelength)**2)
        self.phi_0 = 20*np.sqrt(3)/(self.diameter/self.wavelength)
        
        # Limit parameters
        self.phi_1 = self.phi_0*np.sqrt(20/3)
        self.phi_2 = 10**((49-self.g_0)/25)
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros_like(phi)
        
        # First part
        interval_idx = np.where(np.logical_and(phi >= 0, phi < self.phi_1))
        gain[interval_idx] = self.g_0 - 3*(phi[interval_idx]/self.phi_0)**2
        # Second part
        interval_idx = np.where(np.logical_and(phi >= self.phi_1, phi < self.phi_2))
        gain[interval_idx] = self.g_0 - 20
        # Third part
        interval_idx = np.where(np.logical_and(phi >= self.phi_2, phi < 48))
        gain[interval_idx] = 29 - 25*np.log10(phi[interval_idx])
        # Fourth part
        interval_idx = np.where(np.logical_and(phi >= 48, phi < 80))
        gain[interval_idx] = -13
        # Fifth part
        interval_idx = np.where(np.logical_and(phi >= 80, phi < 120))
        gain[interval_idx] = -8
        # Sixth part
        interval_idx = np.where(np.logical_and(phi >= 120, phi <= 180))
        gain[interval_idx] = -13
        
        return gain
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    par = ParametersRas();
    par.diameter = 1
    par.antenna_efficiency = 1
    par.frequency = 43000
    par.SPEED_OF_LIGHT = 3e8
        
    antenna1 = AntenaSA509(par)
    par.diameter = 7
    antenna7 = AntenaSA509(par)
    par.diameter = 10
    antenna10 = AntenaSA509(par)

    phi = np.linspace(0.1, 180, num = 100000)
    gain1  = antenna1.calculate_gain(phi_vec = phi)
    gain7  = antenna7.calculate_gain(phi_vec = phi)
    gain10  = antenna10.calculate_gain(phi_vec = phi)
    
    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.semilogx(phi, gain1, "-b", label = "$f = 43$ $GHz,$ $D = 1$ $m$")
    plt.semilogx(phi, gain7, "-r", label = "$f = 43$ $GHz,$ $D = 7$ $m$")
    plt.semilogx(phi, gain10, "-k", label = "$f = 43$ $GHz,$ $D = 10$ $m$")
    
    plt.title("ITU-R SA.509-3 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))

    plt.grid()
    plt.show()  
    
    
>>>>>>> dev-ras
