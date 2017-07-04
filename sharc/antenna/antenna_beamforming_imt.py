# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:35:51 2017

@author: Calil
"""

import numpy as np
import matplotlib.pyplot as plt

from sharc.antenna.antenna_element_imt import AntennaElementImt
from sharc.antenna.antenna import Antenna
from sharc.support.named_tuples import AntennaPar
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaBeamformingImt(Antenna):
    """
    Implements an antenna array
    
    Attributes
    ----------
        element (AntennaElementImt): antenna element
        n_rows (int): number of rows in array
        n_cols (int): number of columns in array
        dh (float): horizontal element spacing over wavelenght (d/lambda)
        dv (float): vertical element spacing over wavelenght (d/lambda)
        beams (list): vertical and horizontal tilts of beams
    """
    
    def __init__(self,par: AntennaPar, azimuth: float, elevation: float):
        """
        Constructs an AntennaBeamformingImt object.
        Does not receive angles in local coordinate system.
        Elevation taken with x axis as reference.
        
        Parameters
        ---------
            param (AntennaPar): antenna IMT parameters
            azimuth (float): antenna's physical azimuth inclination
            elevation (float): antenna's physical elevation inclination
                referenced in the x axis
        """
        self.param = par
        
        self.element = AntennaElementImt(par)
        
        self.__azimuth = azimuth
        self.__elevation = elevation
        
        self.__n_rows = par.n_rows
        self.__n_cols = par.n_columns
        self.__dh = par.element_horiz_spacing
        self.__dv = par.element_vert_spacing
        
        self.__beams_list = []
        self.__w_vec_list = []
    
    def add_beam(self, phi_etilt: float, theta_etilt: float):
        """
        Add new beam to antenna.
        Does not receive angles in local coordinate system.
        Theta taken with z axis as reference.
        
        Parameters
        ----------
            phi_etilt (float): azimuth electrical tilt angle [degrees]
            theta_etilt (float): elevation electrical tilt angle [degrees]
        """
        phi, theta = self.to_local_coord(phi_etilt,theta_etilt)
        self.__beams_list.append((phi,theta-90))
        self.__w_vec_list.append(self._weight_vector(phi,theta-90))
        
    def calculate_gain(self,phi_vec: np.array, theta_vec: np.array, beams_l: np.array) -> np.array:
        """
        Calculates the gain in the given direction.
        Does not receive angles in local coordinate system.
        Theta taken with z axis as reference.
        
        Parameters
        ----------
        phi_vec (np.array): azimuth angles [degrees]
        theta_vec (np.array): elevation angles [degrees]
        beam_l (np.array of int): index of beams for gain calculation
            
        Returns
        -------
        gains (np.array): gain corresponding to each of the given directions.
        """
        lo_phi_vec, lo_theta_vec = self.to_local_coord(phi_vec,theta_vec)
        
        n_direct = len(lo_theta_vec)
        
        gains = np.zeros(n_direct)
        
        for g in range(n_direct):
                gains[g] = self._beam_gain(lo_phi_vec[g],lo_theta_vec[g],\
                     beams_l[g])
                
        return gains
    
    def reset_beams(self):
        self.__beams_list = []
        self.__w_vec_list = []
        
    @property
    def azimuth(self):
        return self.__azimuth
    
    @property
    def elevation(self):
        return self.__elevation
    
    @property
    def n_rows(self):
        return self.__n_rows
    
    @property
    def n_cols(self):
        return self.__n_cols
    
    @property
    def dh(self):
        return self.__dh
    
    @property
    def dv(self):
        return self.__dv
    
    @property
    def beams_list(self):
        return self.__beams_list
    
    @property
    def w_vec_list(self):
        return self.__w_vec_list
    
    def _super_position_vector(self,phi: float, theta: float) -> np.array:
        """
        Calculates super position vector.
        Angles are in the local coordinate system.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phi (float): azimuth angle [degrees]
            
        Returns
        -------
            v_vec (np.array): superposition vector
        """
        r_phi = np.deg2rad(phi)
        r_theta = np.deg2rad(theta)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.cos(r_theta) + \
                  (m - 1)*self.dh*np.sin(r_theta)*np.sin(r_phi)
        
        v_vec = np.exp(2*np.pi*1.0j*exp_arg)
        
        return v_vec
        
    def _weight_vector(self, phi_tilt: float, theta_tilt: float) -> np.array:
        """
        Calculates super position vector.
        Angles are in the local coordinate system.
        
        Parameters
        ----------
            phi_tilt (float): electrical horizontal steering [degrees]
            theta_tilt (float): electrical down-tilt steering [degrees]
            
        Returns
        -------
            w_vec (np.array): weighting vector
        """
        r_phi = np.deg2rad(phi_tilt)
        r_theta = np.deg2rad(theta_tilt)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.sin(r_theta) - \
                  (m - 1)*self.dh*np.cos(r_theta)*np.sin(r_phi)
        
        w_vec = (1/np.sqrt(self.n_rows*self.n_cols))*\
                np.exp(2*np.pi*1.0j*exp_arg)
        
        return w_vec        
    
    def _beam_gain(self,phi: float, theta: float, beam = -1) -> float:
        """
        Calculates gain for a single beam in a given direction.
        Angles are in the local coordinate system.
        
        Parameters
        ----------
            phi (float): azimuth angle [degrees]
            theta (float): elevation angle [degrees]
            beam (int): Optional, beam index. If not provided, maximum gain is 
                calculated
            
        Returns
        -------
            gain (float): beam gain [dBi]
        """
        
        element_g = self.element.element_pattern(phi,theta)
        
        v_vec = self._super_position_vector(phi,theta)
        
        if(beam == -1):
            w_vec = self._weight_vector(phi,theta-90)
            array_g = 10*np.log10(abs(np.sum(np.multiply(v_vec,w_vec)))**2)
        else:
            array_g = 10*np.log10(abs(np.sum(np.multiply(v_vec,\
                                            self.__w_vec_list[beam])))**2)
        
        gain = element_g + array_g
        
        return gain
    
    def to_local_coord(self,phi: float, theta: float) -> tuple:
        
        lo_theta = np.ravel(np.array([theta + self.elevation]))
        lo_phi = np.ravel(np.array([phi - self.azimuth]))
        
        ofb_theta = np.where(np.logical_or(lo_theta < 0,lo_theta > 180))
        lo_theta[ofb_theta] = np.mod((360 - lo_theta[ofb_theta]),180)
        lo_phi[ofb_theta] = lo_phi[ofb_theta] + 180
        
        ofb_phi = np.where(np.logical_or(lo_phi < -180,lo_phi > 180))
        lo_phi[ofb_phi] = np.mod(lo_phi[ofb_phi],360)
        ofb_phi = np.where(lo_phi > 180)
        lo_phi[ofb_phi] = lo_phi[ofb_phi] - 360
        
        return lo_phi, lo_theta
    
###############################################################################
class PlotAntennaPattern(object):
    """
    Plots imt antenna pattern.
    """
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
            
        top_y_lim = np.ceil(np.max(gain)/10)*10

        fig = plt.figure(figsize=(15,8))
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
        phi = phi_escan*np.ones(np.size(theta))

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
        if(np.max(gain) > top_y_lim): top_y_lim = np.ceil(np.max(gain)/10)*10
        ax2.set_ylim(top_y_lim - 60,top_y_lim)
        
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
        plt.show()
        
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