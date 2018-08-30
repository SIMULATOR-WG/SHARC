# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 14:12:25 2018

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

class PropagationHDFSSBuildingSide(Propagation):
    """
    
    """
    def __init__(self, param: ParametersFssEs, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        
        self.param = param
        
        # Building dimentions
        self.b_w = 120
        self.b_d = 50
        self.b_h = 3
        self.s_w = 30
        self.b_tol = 0.05
        
        self.LOSS_PER_FLOOR = 50
        
        self.propagation_fspl = PropagationFreeSpace(random_number_gen)
        self.propagation_p1411 = PropagationP1411(random_number_gen,
                                                  above_clutter = False)
        self.building_entry = PropagationBuildingEntryLoss(self.random_number_gen)
        
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        
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
        
        # Define which stations are on the building in front
        next_build = self.is_next_building(imt_x,imt_y,
                                           es_x, es_y)
        not_next_build = np.logical_not(next_build)
        
        # Define which stations are in other buildings
        other_build = np.logical_and(not_same_build,not_next_build)
        
        # Path loss
        loss = np.zeros_like(d)
        
#        # Use a loss per floor
#        loss[:,same_build] += self.get_same_build_loss(imt_z[same_build],
#                                                       es_z)
        loss[:,same_build] += self.propagation_fspl.get_loss(distance_3D=d[:,same_build],
                                                             frequency=f[:,same_build])
        
        loss[:,next_build] += self.propagation_p1411.get_loss(distance_3D=d[:,next_build],
                                                              frequency=f[:,next_build],
                                                              los=True,
                                                              shadow=self.param.shadow_enabled)
        
        loss[:,other_build] += self.propagation_p1411.get_loss(distance_3D=d[:,other_build],
                                                               frequency=f[:,other_build],
                                                               los=False,
                                                               shadow=self.param.shadow_enabled)
    
        # Building entry loss
        if self.param.building_loss_enabled:
            build_loss = self.get_building_loss(imt_sta_type,
                                                f,
                                                elevation)
        else:
            build_loss = 0.0
            
        # Diffraction loss
        diff_loss = np.zeros_like(loss)
                
        # Compute final loss
        loss = loss + build_loss + diff_loss
        
        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)
            
        return loss, build_loss, diff_loss
    
    
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
        building_y_range = (es_y - self.b_d/2) + (1 + self.b_tol)*np.array([-self.b_d/2,+self.b_d/2])
        
        is_in_x = np.logical_and(imt_x >= building_x_range[0],imt_x <= building_x_range[1])
        is_in_y = np.logical_and(imt_y >= building_y_range[0],imt_y <= building_y_range[1])
        
        is_in_building = np.logical_and(is_in_x,is_in_y)
        
        return is_in_building
    
    def get_same_build_loss(self,imt_z,es_z):
        floor_number = imt_z - es_z
        floor_number[floor_number >= 0] = np.floor(floor_number[floor_number >= 0]/self.b_h)
        floor_number[floor_number < 0] = np.ceil(floor_number[floor_number < 0]/self.b_h)
        
        loss = self.LOSS_PER_FLOOR*floor_number
        
        return loss
    
    def is_next_building(self,imt_x,imt_y, es_x, es_y):
        same_building_x_range = es_x + (1 + self.b_tol)*np.array([-self.b_w/2,+self.b_w/2])
        same_building_y_range = (es_y - self.b_d/2) + (1 + self.b_tol)*np.array([-self.b_d/2,+self.b_d/2])
        
        next_building_x_range = same_building_x_range
        next_building_y_range = same_building_y_range + self.b_d + self.s_w
        
        is_in_x = np.logical_and(imt_x >= next_building_x_range[0],imt_x <= next_building_x_range[1])
        is_in_y = np.logical_and(imt_y >= next_building_y_range[0],imt_y <= next_building_y_range[1])
        
        is_in_next_building = np.logical_and(is_in_x,is_in_y)
        
        return is_in_next_building
    
    