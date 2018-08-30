# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 17:28:47 2018

@author: Calil
"""

import numpy as np
import sys
from shapely.geometry import LineString, Polygon, Point

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_p1411 import PropagationP1411
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.support.enumerations import StationType

class PropagationHDFSSRoofTop(Propagation):
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
        
        self.HIGH_LOSS = 4000
        self.LOSS_PER_FLOOR = 50
        
        self.propagation_fspl = PropagationFreeSpace(random_number_gen)
        self.propagation_p1411 = PropagationP1411(random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(self.random_number_gen)
        
        self.SPEED_OF_LIGHT = 299792458.0
        
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
        # Parse entries
        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]
            
        elevation = np.transpose(kwargs["elevation"])
        imt_sta_type = kwargs["imt_sta_type"]
        f = kwargs["frequency"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        
        imt_x = kwargs['imt_x']
        imt_y = kwargs['imt_y']
        imt_z = kwargs['imt_z']
        es_x = kwargs["es_x"]
        es_y = kwargs["es_y"]
        es_z = kwargs["es_z"]
        
        # Define which stations are on the same building
        same_build = self.is_same_building(imt_x,imt_y,
                                           es_x, es_y)
        not_same_build = np.logical_not(same_build)
        
        # Define boolean ranges
        fspl_bool = d <= self.fspl_dist
        
        fspl_to_los_bool = np.logical_and(d > self.fspl_dist,
                                          d <= self.fspl_to_los_dist)
        
        los_bool = np.logical_and(d > self.fspl_to_los_dist,
                                  d <= self.los_dist)
        
        los_to_nlos_bool = np.logical_and(d > self.los_dist,
                                          d <= self.los_to_nlos_dist)
        
        nlos_bool = d > self.los_to_nlos_dist
        
        # Define indexes
        same_build_idx = np.where(same_build)[0]
        fspl_idx = np.where(fspl_bool)[1]
        fspl_to_los_idx = np.where(fspl_to_los_bool)[1]
        los_idx = np.where(los_bool)[1]
        los_to_nlos_idx = np.where(los_to_nlos_bool)[1]
        nlos_idx = np.where(nlos_bool)[1]
        
        # Path loss
        loss = np.zeros_like(d)
        
        if not self.param.same_building_enabled:
            loss[:,same_build_idx] += self.HIGH_LOSS
        else:
            loss[:,same_build_idx] += self.get_same_build_loss(imt_z[same_build_idx],
                                                               es_z)
        
        loss[:,fspl_idx] += self.propagation_fspl.get_loss(distance_3D=d[:,fspl_idx],
                                                           frequency=f[:,fspl_idx])
        loss[:,fspl_to_los_idx] += self.interpolate_fspl_to_los(d[:,fspl_to_los_idx],
                                                                f[:,fspl_to_los_idx],
                                                                self.param.shadow_enabled)
        loss[:,los_idx] += self.propagation_p1411.get_loss(distance_3D=d[:,los_idx],
                                                           frequency=f[:,los_idx],
                                                           los=True,
                                                           shadow=self.param.shadow_enabled)
        loss[:,los_to_nlos_idx] += self.interpolate_los_to_nlos(d[:,los_to_nlos_idx],
                                                                f[:,los_to_nlos_idx],
                                                                self.param.shadow_enabled)
        loss[:,nlos_idx] += self.propagation_p1411.get_loss(distance_3D=d[:,nlos_idx],
                                                            frequency=f[:,nlos_idx],
                                                            los=False,
                                                            shadow=self.param.shadow_enabled)
    
        # Building entry loss
        if self.param.building_loss_enabled:
            build_loss = np.zeros_like(loss)
            build_loss[0,not_same_build] = self.get_building_loss(imt_sta_type,
                                                                  f[:,not_same_build],
                                                                  elevation[:,not_same_build])
        else:
            build_loss = 0.0
            
        # Diffraction loss
        diff_loss = np.zeros_like(loss)
        if self.param.diffraction_enabled:
            h, d1, d2 = self.get_diff_distances(imt_x,imt_y, imt_z, 
                                                 es_x, es_y,  es_z)
            diff_loss = np.zeros_like(loss)
            diff_loss[0,not_same_build] = self.get_diffraction_loss(h[not_same_build],
                                                                    d1[not_same_build],
                                                                    d2[not_same_build],
                                                                    f[:,not_same_build])
                
        # Compute final loss
        loss = loss + build_loss + diff_loss
        
        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)
            
        return loss, build_loss, diff_loss
    
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
        
        building_x_range = es_x + (1 + self.b_tol)*np.array([-self.b_w/2,+self.b_w/2])
        building_y_range = es_y + (1 + self.b_tol)*np.array([-self.b_d/2,+self.b_d/2])
        
        is_in_x = np.logical_and(imt_x >= building_x_range[0],imt_x <= building_x_range[1])
        is_in_y = np.logical_and(imt_y >= building_y_range[0],imt_y <= building_y_range[1])
        
        is_in_building = np.logical_and(is_in_x,is_in_y)
        
        return is_in_building
    
    def get_same_build_loss(self,imt_z,es_z):
        floor_number = np.floor_divide((es_z - imt_z),3) + 1
        
        loss = self.LOSS_PER_FLOOR*floor_number
        
        return loss
    
    def get_diff_distances(self,imt_x,imt_y, imt_z, es_x, es_y, es_z, dist_2D=False):
        
        build_poly = Polygon([[es_x + self.b_w/2, es_y + self.b_d/2],
                              [es_x - self.b_w/2, es_y + self.b_d/2],
                              [es_x - self.b_w/2, es_y - self.b_d/2],
                              [es_x + self.b_w/2, es_y - self.b_d/2]])
        es_point = Point([es_x,es_y])
        
        d1_2D = np.zeros_like(imt_x)
        d2_2D = np.zeros_like(imt_x)
        dist = np.zeros_like(imt_x)
        for k,(x,y) in enumerate(zip(imt_x,imt_y)):
            line = LineString([[es_x,es_y],[x,y]])
            intersection_line = line.intersection(build_poly)
            d1_2D[k] = intersection_line.length
            
            imt_point = Point([x,y])
            dist[k] = es_point.distance(imt_point)
            d2_2D[k] = dist[k] - d1_2D[k]
            
        if dist_2D:
            return d1_2D, d2_2D
        
        build_height = es_z - 1
        z_in_build = imt_z + (es_z - imt_z)*d2_2D/dist
        h = build_height - z_in_build
        
        d1 = np.sqrt(1**2 + np.power(d1_2D,2))
        d2 = np.sqrt(np.power((build_height - imt_z),2) + np.power(d2_2D,2))
        
        return h, d1, d2
    
    def get_diffraction_loss(self,h, d1, d2, f):
        
        wavelength =  self.SPEED_OF_LIGHT/(f*1e6)
        
        v = h*np.sqrt((2/wavelength)*(1/d1 + 1/d2))
        
        loss = np.zeros_like(v)
        
        v_idx = v > -0.7
        
        loss[v_idx] = 6.9 + 20*np.log10(np.sqrt(np.power((v[v_idx] - 0.1), 2) + 1)\
                                        + v[v_idx] - 0.1) 
        return loss
    
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
    d = np.array([list(d)])
    f = 40e3*np.ones_like(d)
    ele = np.transpose(np.zeros_like(d))
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
    
    ravel_d = np.ravel(d)
    ravel_loss_interp = np.ravel(loss_interp)
    ravel_loss_no_interp = np.ravel(loss_no_interp)
    plt.plot(ravel_d,ravel_loss_interp,'k-',label='Interpolated')
    plt.plot(ravel_d,ravel_loss_no_interp,'k--',label='Not Interpolated')
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
    
    ravel_loss = np.ravel(loss)
    plt.plot(ravel_d,ravel_loss)
    plt.xlabel("Distance [m]")
    plt.ylabel("Path Loss [dB]")
    plt.grid()
    plt.show()
    
    