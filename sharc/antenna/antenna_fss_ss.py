# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:16:29 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import sys

class AntennaFssSs(Antenna):
    """
    Implements the antenna pattern for FSS space station according to Report on 
    the second meeting of Task Group 5/1 (Geneva, Switzerland, 15-23 May 2017).
    
    Gain to use in the sidelobe (for angles beyond aψ0 °): Same as Recommendation 
    ITU-R S.672-4 with LS=-25 dB or with the value provided by WP 4A and a 
    reduction of 3dB.
    
    Curve to use inside the main lobe (i.e. from 0° to aψ0 °): Curve according 
    to pattern APSREC408V01 from the BR space station antenna pattern library. 
    
    """
    
    def __init__(self, param: ParametersFssSs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.l_s = param.antenna_l_s
        
        if self.l_s == -20:
            self.a = 2.58
        elif self.l_s == -25:
            self.a = 2.88
        elif self.l_s == -30:
            self.a = 3.16
        else:
            sys.stderr.write("ERROR\nInvalid AntennaFssSs L_s parameter: " + self.l_s)
            sys.exit(1)

        self.b = 6.32
        
        self.psi_0 = param.antenna_3_dB/2
        self.psi_1 = self.psi_0 * np.power(10, (self.peak_gain + self.l_s + 20)/25)
            
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros(len(psi))
        
        idx_1 = np.where(psi <= self.a * self.psi_0)[0]
        gain[idx_1] = self.peak_gain - 12 * np.power(psi[idx_1]/(2*self.psi_0), 2)
            
        idx_2 = np.where((self.a * self.psi_0 < psi) & (psi <= self.b * self.psi_0 ))[0]
        gain[idx_2] = self.peak_gain + self.l_s - 3
        
        idx_3 = np.where((self.b * self.psi_0 < psi) & (psi <= self.psi_1))[0]
        gain[idx_3] = self.peak_gain + self.l_s + 20 - 25 * np.log10(psi[idx_3]/self.psi_0) - 3
        
        return gain
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from sharc.antenna.antenna_s672 import AntennaS672

    # initialize antenna parameters
    param = ParametersFssSs()
    param.antenna_gain = 50
    param.antenna_pattern = "FSS_SS"
    param.antenna_3_dB = 2
    psi = np.linspace(0.1, 30, num = 1000)
    
    param.antenna_l_s = -25    
    antenna_fss = AntennaFssSs(param)    
    gain_fss = antenna_fss.calculate_gain(phi_vec=psi)

    antenna_672 = AntennaS672(param)    
    gain_672 = antenna_672.calculate_gain(phi_vec=psi)
    
    fig = plt.figure(figsize=(12,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.semilogx(psi, gain_672 - param.antenna_gain, "-b", label="ITU-R S.672-4")
    plt.semilogx(psi, gain_fss - param.antenna_gain, "-r", label="FSS SS")

    plt.ylim((-40, 0.3))
    plt.xlim((0.1, 100))    
    plt.title("FSS space station antenna pattern")
    plt.xlabel("Relative off-axis angle, $\psi/\psi_0$")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="upper right")

    ax = plt.gca()
    ax.set_yticks([-40, -30, -20, -10, 0])
    ax.set_xticks(np.linspace(0.1, 0.9, 9).tolist() + np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()        