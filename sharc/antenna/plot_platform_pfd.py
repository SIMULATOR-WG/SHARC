# -*- coding: utf-8 -*-
"""
Created on Mon May  7 22:43:22 2018

@author: Edgar Souza
"""

import numpy as np
import matplotlib.pyplot as plt

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.support.named_tuples import AntennaPar
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

if __name__ == '__main__':

    param = ParametersAntennaImt()
    param.normalization = False    
    param.bs_normalization_file = "beamforming_normalization\\norm_bs_8x8.npz"    
    param.ue_normalization_file = "beamforming_normalization\\norm_ue_4x4.npz"    

    param.bs_element_pattern = "M2101"
    param.bs_tx_element_max_g    = 5
    param.bs_tx_element_phi_deg_3db  = 65
    param.bs_tx_element_theta_deg_3db = 65
    param.bs_tx_element_am       = 30
    param.bs_tx_element_sla_v    = 30
    param.bs_tx_n_rows           = 8
    param.bs_tx_n_columns        = 8
    param.bs_tx_element_horiz_spacing = 0.5
    param.bs_tx_element_vert_spacing = 0.5
    param.bs_downtilt_deg = 0

    param.ue_element_pattern = "M2101"
    param.ue_tx_element_max_g    = 5
    param.ue_tx_element_phi_deg_3db  = 90
    param.ue_tx_element_theta_deg_3db = 90
    param.ue_tx_element_am       = 25
    param.ue_tx_element_sla_v    = 25
    param.ue_tx_n_rows           = 4
    param.ue_tx_n_columns        = 4
    param.ue_tx_element_horiz_spacing = 0.5
    param.ue_tx_element_vert_spacing = 0.5

    par_bs = param.get_antenna_parameters("BS", "TX")
    ant_bs = AntennaBeamformingImt(par_bs, 0, 0)
    par_ue = param.get_antenna_parameters("UE", "TX")
    ant_ue = AntennaBeamformingImt(par_ue, 0, 0)
    
    param.normalization = True
    par_bs_norm = param.get_antenna_parameters("BS", "TX")
    ant_bs_norm = AntennaBeamformingImt(par_bs_norm, 0, 0)
    par_ue_norm = param.get_antenna_parameters("UE", "TX")
    ant_ue_norm = AntennaBeamformingImt(par_ue_norm, 0, 0)
    
    phi_escan = 0
    theta_tilt = 90
    ant_bs.add_beam(phi_escan, theta_tilt)
    ant_bs_norm.add_beam(phi_escan, theta_tilt)
    ant_ue.add_beam(phi_escan, theta_tilt)
    ant_ue_norm.add_beam(phi_escan, theta_tilt)

    theta = np.linspace(0, 180, num = 1000)
    phi = phi_escan*np.ones(np.size(theta))

    gain_bs = ant_bs.calculate_gain(phi_vec = phi,
                                    theta_vec = theta,
                                    beams_l = np.zeros_like(phi, dtype=int))
    gain_bs_norm = ant_bs_norm.calculate_gain(phi_vec = phi,
                                              theta_vec = theta,
                                              beams_l = np.zeros_like(phi, dtype=int))
    gain_ue = ant_ue.calculate_gain(phi_vec = phi,
                                    theta_vec = theta,
                                    beams_l = np.zeros_like(phi, dtype=int))
    gain_ue_norm = ant_ue_norm.calculate_gain(phi_vec = phi,
                                              theta_vec = theta,
                                              beams_l = np.zeros_like(phi, dtype=int))

    protection_criteria = -6
    frequency = 27250*1e6
    wavelength = 299792458/frequency
    thermal_noise = -144
    noise_figure = 10
    
    pfd_bs = protection_criteria \
                + 10*np.log10(4*np.pi/(wavelength**2)) \
                + thermal_noise \
                - gain_bs \
                + noise_figure
    pfd_bs_norm = protection_criteria \
                    + 10*np.log10(4*np.pi/(wavelength**2)) \
                    + thermal_noise \
                    - gain_bs_norm \
                    + noise_figure
    pfd_ue = protection_criteria \
                + 10*np.log10(4*np.pi/(wavelength**2)) \
                + thermal_noise \
                - gain_ue \
                + noise_figure
    pfd_ue_norm = protection_criteria \
                    + 10*np.log10(4*np.pi/(wavelength**2)) \
                    + thermal_noise \
                    - gain_ue_norm \
                    + noise_figure
                
    plt.figure(figsize = (15, 5), facecolor = "w", edgecolor = "k")
    ax = plt.subplot(121)
    
    ax.plot(theta, pfd_bs, "-b", label = "BS", linewidth = 1)
    ax.plot(theta, pfd_bs_norm, "--b", label = "BS (norm. antenna)", linewidth = 1)
    ax.plot(theta, pfd_ue, "-g", label = "UE", linewidth = 1)
    ax.plot(theta, pfd_ue_norm, "--g", label = "UE (norm. antenna)", linewidth = 1)
        
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel("PFD mask [$dB(W/m^2)$ in 1 MHz]")
    ax.set_xlim((0, 180))
    ax.set_ylim((-120, -60))
    ax.legend(bbox_to_anchor = (1.03, 1), loc = 2, borderaxespad = 0., fontsize = 12)
    ax.grid(True)    
    plt.show()