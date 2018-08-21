# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 17:28:47 2018

@author: Calil
"""

import numpy as np
import sys

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_p1411 import PropagationP1411
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.support.enumerations import StationType

class PropagationHDFSS(Propagation):
    """
    This is a wrapper class which can be used for indoor simulations. It
    calculates the basic path loss between IMT stations and HDFSS Earth Stations.
    It uses:
        FSPL for distance < 55m
        P.1411 LOS for distances 55m < distance < 260m
        P.1411 NLOS for distances distance > 260m
    """
    def __init__(self, param: ParametersFssEs, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        
        self.param = param
        
        self.fspl_dist = 35
        self.fspl_to_los_dist = 55
        self.los_dist = 100
        self.los_to_nlos_dist = 260
        
        # Building dimentions
        self.b_w = 120
        self.b_d = 50
        self.b_tol = 0.05
        
        self.HIGH_LOSS = 400
        
        self.propagation_fspl = PropagationFreeSpace(random_number_gen)
        self.propagation_p1411 = PropagationP1411(random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(self.random_number_gen)
        
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for given distances and frequencies

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequencie [MHz]
            indoor (np.array) : indicates whether UE is indoor
            shadowing (bool) : if shadowing should be added or not
            number_of_sectors (int): number of sectors in a node (default = 1)

        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]
            
        elevation = np.transpose(kwargs["elevation"])
        imt_sta_type = kwargs["imt_sta_type"]
        f = kwargs["frequency"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        
        if not self.param.same_building_enabled:
            imt_x = kwargs['imt_x']
            imt_y = kwargs['imt_y']
            imt_z = kwargs['imt_z']
            es_x = kwargs["es_x"]
            es_y = kwargs["es_y"]
            es_z = kwargs["es_z"]
            same_build = self.is_same_building(imt_x,imt_y,
                                               es_x, es_y)
        else:
            same_build = np.zeros_like(d, dtype=bool)  
        not_same_build = np.logical_not(same_build)
        
        fspl_bool = np.logical_and(d <= self.fspl_dist,not_same_build)
        
        fspl_to_los_bool = np.logical_and(d > self.fspl_dist,
                                          d <= self.fspl_to_los_dist)
        fspl_to_los_bool = np.logical_and(fspl_to_los_bool,not_same_build)
        
        los_bool = np.logical_and(d > self.fspl_to_los_dist,d <= self.los_dist)
        los_bool = np.logical_and(los_bool,not_same_build)
        
        los_to_nlos_bool = np.logical_and(d > self.los_dist,d <= self.los_to_nlos_dist)
        los_to_nlos_bool = np.logical_and(los_to_nlos_bool,not_same_build)
        
        nlos_bool = np.logical_and(d > self.los_to_nlos_dist,not_same_build)
        
        fspl_idx = np.where(fspl_bool)
        fspl_to_los_idx = np.where(fspl_to_los_bool)
        los_idx = np.where(los_bool)
        los_to_nlos_idx = np.where(los_to_nlos_bool)
        nlos_idx = np.where(nlos_bool)
        
        loss = np.zeros_like(d)
        
        loss[fspl_idx] = self.propagation_fspl.get_loss(distance_3D=d[fspl_idx],
                                                        frequency=f[fspl_idx])
        loss[fspl_to_los_idx] = self.interpolate_fspl_to_los(d[fspl_to_los_idx],
                                                             f[fspl_to_los_idx],
                                                             self.param.shadow_enabled)
        loss[los_idx] = self.propagation_p1411.get_loss(distance_3D=d[los_idx],
                                                        frequency=f[los_idx],
                                                        los=True,
                                                        shadow=self.param.shadow_enabled)
        loss[los_to_nlos_idx] = self.interpolate_los_to_nlos(d[los_to_nlos_idx],
                                                             f[los_to_nlos_idx],
                                                             self.param.shadow_enabled)
        loss[nlos_idx] = self.propagation_p1411.get_loss(distance_3D=d[nlos_idx],
                                                         frequency=f[nlos_idx],
                                                         los=False,
                                                         shadow=self.param.shadow_enabled)
    
        if self.param.building_loss_enabled:
            build_loss = self.get_building_loss(imt_sta_type,f,elevation)
        else:
            build_loss = 0.0
                
        loss = loss + build_loss
        
        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)
            
        return loss
    
    def interpolate_fspl_to_los(self,dist,freq,shad): 
        fspl_loss = self.propagation_fspl.get_loss(distance_3D=self.fspl_dist,
                                                   frequency=freq)
        los_loss = self.propagation_p1411.get_loss(distance_3D=self.fspl_to_los_dist,
                                                   frequency=freq,
                                                   los=True,
                                                   shadow=False)
        
        loss = (dist - self.fspl_dist)*(los_loss - fspl_loss)/(self.fspl_to_los_dist - self.fspl_dist) + fspl_loss
        
        if shad:
            interp_sigma = (dist - self.fspl_dist)*(self.propagation_p1411.los_sigma)/(self.fspl_to_los_dist - self.fspl_dist)
            loss = loss + self.random_number_gen.normal(0.0,interp_sigma)
            
        return loss
    
    def interpolate_los_to_nlos(self,dist,freq,shad): 
        los_loss = self.propagation_p1411.get_loss(distance_3D=self.los_dist,
                                                   frequency=freq,
                                                   los=True,
                                                   shadow=False)
        nlos_loss = self.propagation_p1411.get_loss(distance_3D=self.los_to_nlos_dist,
                                                    frequency=freq,
                                                    los=False,
                                                    shadow=False)
        
        loss = (dist-self.los_dist)*(nlos_loss - los_loss)/(self.los_to_nlos_dist - self.los_dist) + los_loss
        
        if shad:
            interp_sigma = (dist-self.los_dist)*(self.propagation_p1411.nlos_sigma - self.propagation_p1411.los_sigma)/(self.los_to_nlos_dist - self.los_dist) +\
                           self.propagation_p1411.los_sigma
            loss = loss + self.random_number_gen.normal(0.0,interp_sigma)
            
        return loss
    
    def get_building_loss(self,imt_sta_type,f,elevation):
        if imt_sta_type is StationType.IMT_UE:
            build_loss = self.building_entry.get_loss(f, elevation)
        elif imt_sta_type is StationType.IMT_BS:
            if self.param.bs_building_entry_loss_type == 'P2109_RANDOM':
                build_loss = self.building_entry.get_loss(f, elevation)
            elif self.param.bs_building_entry_loss_type == 'P2109_FIXED':
                build_loss = self.building_entry.get_loss(f, elevation, prob=self.param.bs_building_entry_loss_prob)
            elif self.param.bs_building_entry_loss_type == 'FIXED_VALUE':
                build_loss = self.param.bs_building_entry_loss_value
            else:
                sys.stderr.write("ERROR\nBuilding entry loss type: " + 
                                 self.param.bs_building_entry_loss_type)
                sys.exit(1)
                
        return build_loss
    
    def is_same_building(self,imt_x,imt_y, es_x, es_y):
        
        building_x_range = (1 + self.b_tol)*np.array([es_x - self.w/2, es_x + self.w/2])
        building_y_range = (1 + self.b_tol)*np.array([es_y - self.d/2, es_y + self.d/2])
        
        is_in_x = np.logical_and(imt_x >= building_x_range[0],imt_x <= building_x_range[1])
        is_in_y = np.logical_and(imt_y >= building_y_range[0],imt_y <= building_y_range[1])
        
        is_in_building = np.logical_and(is_in_x,is_in_y)
        
        return is_in_building
    
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    rnd = np.random.RandomState(101)
    par = ParametersFssEs()
    par.building_loss_enabled = False
    par.shadow_enabled = False
    par.same_building_enabled = True
    par.bs_building_entry_loss_type = 'FIXED_VALUE'
    par.bs_building_entry_loss_prob = 0.5
    par.bs_building_entry_loss_value = 50
    prop = PropagationHDFSS(par,rnd)
    
    d = np.linspace(5,1000,num=2000)
    f = 40e3*np.ones_like(d)
    ele = np.zeros_like(d)
    sta_type = StationType.IMT_BS
        
    # Without shadow
    loss_interp = prop.get_loss(distance_3D=d,
                         frequency=f,
                         elevation=ele,
                         imt_sta_type=sta_type)
    prop.fspl_dist = prop.fspl_to_los_dist
    prop.los_dist = prop.los_to_nlos_dist
    loss_no_interp = prop.get_loss(distance_3D=d,
                                   frequency=f,
                                   elevation=ele,
                                   imt_sta_type=sta_type)
    
    plt.plot(d,loss_interp,'k-',label='Interpolated')
    plt.plot(d,loss_no_interp,'k--',label='Not Interpolated')
    plt.legend()
    plt.xlabel("Distance [m]")
    plt.ylabel("Path Loss [dB]")
    plt.grid()
    plt.show()
    
    # With shadow
    par.shadow_enabled = True
    prop = PropagationHDFSS(par,rnd)
    loss = prop.get_loss(distance_3D=d,
                         frequency=f,
                         elevation=ele,
                         imt_sta_type=sta_type)
    
    plt.plot(d,loss)
    plt.xlabel("Distance [m]")
    plt.ylabel("Path Loss [dB]")
    plt.grid()
    plt.show()